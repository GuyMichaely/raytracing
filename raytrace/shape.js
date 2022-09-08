/* Intersection structure:
 * t:        ray parameter (float), i.e. distance of intersection point to ray's origin
 * position: position (THREE.Vector3) of intersection point
 * normal:   normal (THREE.Vector3) of intersection point
 * material: material of the intersection object
 */
class Intersection {
	constructor() {
		this.t = 0;
		this.position = new THREE.Vector3();
		this.normal = new THREE.Vector3();
		this.material = null;
	}
	set(isect) {
		this.t = isect.t;
		this.position = isect.position;
		this.normal = isect.normal;
		this.material = isect.material;
	}
}

/* Plane shape
 * P0: a point (THREE.Vector3) that the plane passes through
 * n:  plane's normal (THREE.Vector3)
 */
class Plane {
	constructor(P0, n, material) {
		this.P0 = P0.clone();
		this.n = n.clone();
		this.n.normalize();
		this.material = material;
	}
	// Given ray and range [tmin,tmax], return intersection point.
	// Return null if no intersection.
	intersect(ray, tmin, tmax) {
		let temp = this.P0.clone();
		temp.sub(ray.o); // (P0-O)
		let denom = ray.d.dot(this.n); // d.n
		if(denom==0) { return null;	}
		let t = temp.dot(this.n)/denom; // (P0-O).n / d.n
		if(t<tmin || t>tmax) return null; // check range
		let isect = new Intersection();   // create intersection structure
		isect.t = t;
		isect.position = ray.pointAt(t);
		isect.normal = this.n;
		isect.material = this.material;
		return isect;
	}
}

/* Sphere shape
 * C: center of sphere (type THREE.Vector3)
 * r: radius
 */
class Sphere {
	constructor(C, r, material) {
		this.C = C.clone();
		this.r = r;
		this.r2 = r*r;
		this.material = material;
	}
	intersect(ray, tmin, tmax) {
// ===YOUR CODE STARTS HERE===
		const B = 2 * ray.o.clone().sub(this.C).dot(ray.d)
		const C = Math.abs(ray.o.clone().sub(this.C).lengthSq() - this.r2)
		const discriminant_squared = B * B - 4 * C
		if (discriminant_squared < 0) {
			return null;
		}
		const solutions = [-1, 1].map(coef => (-B + coef * Math.pow(discriminant_squared, .5)) / 2).filter(t => t > tmin && t < tmax)

		if (solutions.length) {
			let isect = new Intersection();
			isect.t = Math.min(...solutions);
			isect.position = ray.pointAt(isect.t);
			isect.normal = isect.position.clone().sub(this.C).normalize();
			isect.material = this.material;
			return isect;
		}
		return null

// ---YOUR CODE ENDS HERE---
	}
}

function rrefMatrix(mat) {
    let lead = 0;
    for (let r = 0; r < 3; r++) {
        if (4 <= lead) {
            return;
        }
        let i = r;
        while (mat[i][lead] == 0) {
            i++;
            if (3 == i) {
                i = r;
                lead++;
                if (4 == lead) {
                    return;
                }
            }
        }
 
        let tmp = mat[i];
        mat[i] = mat[r];
        mat[r] = tmp;
 
        let val = mat[r][lead];
        for (let j = 0; j < 4; j++) {
            mat[r][j] /= val;
        }
 
        for (let i = 0; i < 3; i++) {
            if (i == r) continue;
            val = mat[i][lead];
            for (let j = 0; j < 4; j++) {
                mat[i][j] -= val * mat[r][j];
            }
        }
        lead++;
    }
    return mat;
}

function transposeMatrix(mat) {
	let result = [[], [], []];
	for (let i = 0; i < 4; i++) {
		for (let j = 0; j < 3; j++) {
			result[j][i] = mat[i][j];
		}
	}
	return result;
}

class Triangle {
	/* P0, P1, P2: three vertices (type THREE.Vector3) that define the triangle
	 * n0, n1, n2: normal (type THREE.Vector3) of each vertex */
	constructor(P0, P1, P2, material, n0, n1, n2) {
		this.P0 = P0.clone();
		this.P1 = P1.clone();
		this.P2 = P2.clone();
		this.material = material;
		if(n0) this.n0 = n0.clone();
		if(n1) this.n1 = n1.clone();
		if(n2) this.n2 = n2.clone();

		// below you may pre-compute any variables that are needed for intersect function
		// such as the triangle normal etc.
// ===YOUR CODE STARTS HERE===

// ---YOUR CODE ENDS HERE---
	} 

	intersect(ray, tmin, tmax) {
// ===YOUR CODE STARTS HERE===
		// construct linear system and solve it
		const p2minusp0 = this.P2.clone().sub(this.P0)
		const p2minusp1 = this.P2.clone().sub(this.P1)
		const p2minuso 	= this.P2.clone().sub(ray.o)
		const colOrderMatrix = [
			[ray.direction().x, ray.direction().y, ray.direction().z],
			[p2minusp0.x, p2minusp0.y, p2minusp0.z],
			[p2minusp1.x, p2minusp1.y, p2minusp1.z],
			[p2minuso.x, p2minuso.y, p2minuso.z]
		]
		const rowOrderMatrix = transposeMatrix(colOrderMatrix)
		const rrefMat = rrefMatrix(rowOrderMatrix)

		// check that system was solvable
		for (let i = 0; i < 3; i++) {
			if (!rrefMat[i].slice(0, 3).every((val, index) => index == i ? val == 1 : val == 0)) {
				return null;
			}
		}

		// check conditions for intersection
		const solution = rrefMat.map(r => r[3])
		if (solution.slice(1).some(val => val < 0) || solution[1] + solution[2] > 1 || solution[0] < tmin || solution[0] > tmax) {
			return null
		}

		// build Intersection object
		let isect = new Intersection();
		isect.t = Math.min(solution[0]);
		isect.position = ray.pointAt(isect.t);
		isect.material = this.material;
		// compute normal
		if (this.n0 && this.n1 && this.n2) {
			isect.normal = this.n0.clone().multiplyScalar(solution[1]).add(this.n1.clone().multiplyScalar(solution[2])).add(this.n2.clone().multiplyScalar(1 - solution[1] - solution[2])).normalize();
		}
		else {
		const v1 = this.P1.clone().sub(this.P0)
		const v2 = this.P2.clone().sub(this.P1)
		v1.cross(v2)
		isect.normal = v1.normalize()
		}
		
		return isect;

// ---YOUR CODE ENDS HERE---
	}
}

function shapeLoadOBJ(objstring, material, smoothnormal) {
	loadOBJFromString(objstring, function(mesh) { // callback function for non-blocking load
		if(smoothnormal) mesh.computeVertexNormals();
		for(let i=0;i<mesh.faces.length;i++) {
			let p0 = mesh.vertices[mesh.faces[i].a];
			let p1 = mesh.vertices[mesh.faces[i].b];
			let p2 = mesh.vertices[mesh.faces[i].c];
			if(smoothnormal) {
				let n0 = mesh.faces[i].vertexNormals[0];
				let n1 = mesh.faces[i].vertexNormals[1];
				let n2 = mesh.faces[i].vertexNormals[2];
				shapes.push(new Triangle(p0, p1, p2, material, n0, n1, n2));
			} else {
				shapes.push(new Triangle(p0, p1, p2, material));
			}
		}
	}, function() {}, function() {});
}

/* ========================================
 * You can define additional Shape classes,
 * as long as each implements intersect function.
 * ======================================== */
