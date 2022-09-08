/* Ray class:
 * o: origin (THREE.Vector3)
 * d: normalized direction (THREE.Vector3)
 */
class Ray {
	constructor(origin, direction) {
		this.o = origin.clone();
		this.d = direction.clone();
		this.d.normalize();
	}
	pointAt(t) {
		// P(t) = o + t*d
		let point = this.o.clone();
		point.addScaledVector(this.d, t);
		return point;
	}
	direction() { return this.d; }
	origin() { return this.o; }
}

function render() {
	// create canvas of size imageWidth x imageHeight and add to DOM
	let canvas = document.createElement('canvas');
	canvas.width = imageWidth;
	canvas.height = imageHeight;
	canvas.style = 'background-color:red';
	document.body.appendChild(canvas);
	let ctx2d = canvas.getContext('2d'); // get 2d context
	let image = ctx2d.getImageData(0, 0, imageWidth, imageHeight); // get image data
	let pixels = image.data; // get pixel array

	let row=0;
	let idx=0;
	let chunksize=10; // render 10 rows at a time
	console.log('Raytracing started...');
	(function chunk() {
		// render a chunk of rows
		for(let j=row;j<row+chunksize && j<imageHeight;j++) {
			for(let i=0;i<imageWidth;i++,idx+=4) { // i loop
				// compute normalized pixel coordinate (x,y)
				let x = i/imageWidth;
				let y = (imageHeight-1-j)/imageHeight;
				let ray = camera.getCameraRay(x,y);
				let color = raytracing(ray, 0);
				setPixelColor(pixels, idx, color);
			}
		}
		row+=chunksize;  // non-blocking j loop
		if(row<imageHeight) {
			setTimeout(chunk, 0);
			ctx2d.putImageData(image, 0, 0); // display intermediate image
		} else {
			ctx2d.putImageData(image, 0, 0); // display final image
			console.log('Done.')
		}
	})();
}

/* Trace ray in the scene and return color of ray. 'depth' is the current recursion depth.
 * If intersection material has non-null kr or kt, perform recursive ray tracing. */
function raytracing(ray, depth) {
// ===YOUR CODE STARTS HERE===
	let intersect = rayIntersectScene(ray);
	if (intersect == null) {
		return backgroundColor
	}
	if (depth > maxDepth || [intersect.material.kr, intersect.material.kt].every(v => v == null)) {
		return shading(ray, intersect);
	}
	
	let color = new THREE.Color(0, 0, 0);
	if (intersect.material.kr) {
		let reflectedRay = new Ray(intersect.position, reflect(ray.d.clone().multiplyScalar(-1), intersect.normal));
		let reflectedColor = raytracing(reflectedRay, depth + 1);
		color.add(reflectedColor.multiply(intersect.material.kr).add(shading(ray, intersect)))
	}
	if (intersect.material.kt) {
		let refractedRay = new Ray(intersect.position, refract(ray.d, intersect.normal, intersect.material.ior));
		let refractedColor = raytracing(refractedRay, depth + 1);
		color.add(refractedColor.multiply(intersect.material.kt).add(shading(ray, intersect)))
	}
	return color

// ---YOUR CODE ENDS HERE---
}

function getDiffuseReducer(isect) {
	return function diffuseReducer(acc, lightSample) {
		return acc.add(
			lightSample.intensity.clone()
			.multiply(isect.material.kd)
			.multiplyScalar(Math.max(0, isect.normal.dot(lightSample.direction)))
		)
	}
}

function getSpecularReducer(isect, ray) {
	return function specularReducer(acc, lightSample) {
		const reflected = reflect(lightSample.direction, isect.normal)
		return acc.add(
			lightSample.intensity.clone()
			.multiply(isect.material.ks)
			.multiplyScalar(Math.pow(
				Math.max(0, -reflected.dot(ray.d)),
				isect.material.p
			))
		)
	}
}

/* Compute and return shading color given a ray and the intersection point structure. */
function shading(ray, isect) {
	let color = new THREE.Color(0,0,0);
// ===YOUR CODE STARTS HERE===
	// calculate color from ambient light
	if (isect.material.ka != null) {
		color.add(ambientLight.clone().multiply(isect.material.ka));
	}

	const lightSamples = lights.map(light => light.getLight(isect.position)).filter(lightSample => {
		const shadowRay = new Ray(isect.position, lightSample.direction)
		const shadowIsect = rayIntersectScene(shadowRay)
		if (shadowIsect == null) {
			return true
		}
		const distToLight = lightSample.position.clone().sub(isect.position).length()
		const distToIntersection = shadowIsect.position.clone().sub(isect.position).length()
		return distToLight < distToIntersection
	})

	// diffuse
	if (isect.material.kd != null) {
		color.add(lightSamples.reduce(getDiffuseReducer(isect), new THREE.Color(0,0,0)))
	}

	// specular
	if (isect.material.ks != null) {
		color.add(lightSamples.reduce(getSpecularReducer(isect, ray), new THREE.Color(0,0,0)))
	}
// ---YOUR CODE ENDS HERE---
	return color;
}

/* Compute intersection of ray with scene shapes.
 * Return intersection structure (null if no intersection). */
function rayIntersectScene(ray) {
	let tmax = Number.MAX_VALUE;
	let isect = null;
	for(let i=0;i<shapes.length;i++) {
		let hit = shapes[i].intersect(ray, 0.0001, tmax);
		if(hit != null) {
			tmax = hit.t;
			if(isect == null) isect = hit; // if this is the first time intersection is found
			else isect.set(hit); // update intersection point
		}
	}
	return isect;
}

/* Compute reflected vector, by mirroring l around n. */
function reflect(l, n) {
	// r = 2(n.l)*n-l
	let r = n.clone();
	r.multiplyScalar(2*n.dot(l));
	r.sub(l);
	return r;
}

/* Compute refracted vector, given l, n and index_of_refraction. */
function refract(l, n, ior) {
	let mu = (n.dot(l) < 0) ? 1/ior:ior;
	let cosI = l.dot(n);
	let sinI2 = 1 - cosI*cosI;
	if(mu*mu*sinI2>1) return null;
	let sinR = mu*Math.sqrt(sinI2);
	let cosR = Math.sqrt(1-sinR*sinR);
	let r = n.clone();
	if(cosI > 0) {
		r.multiplyScalar(-mu*cosI+cosR);
		r.addScaledVector(l, mu);
	} else {
		r.multiplyScalar(-mu*cosI-cosR);
		r.addScaledVector(l, mu);
	}
	r.normalize();
	return r;
}

/* Convert floating-point color to integer color and assign it to the pixel array. */
function setPixelColor(pixels, index, color) {
	pixels[index+0]=pixelProcess(color.r);
	pixels[index+1]=pixelProcess(color.g);
	pixels[index+2]=pixelProcess(color.b);
	pixels[index+3]=255; // alpha channel is always 255*/
}

/* Multiply exposure, clamp pixel value, then apply gamma correction. */
function pixelProcess(value) {
	value*=exposure; // apply exposure
	value=(value>1)?1:value;
	value = Math.pow(value, 1/2.2);	// 2.2 gamma correction
	return value*255;
}
