#include "flyscene.hpp"
#include <GLFW/glfw3.h>

void Flyscene::initialize(int width, int height) {
  std::cout << "load obj\n";
  // initiliaze the Phong Shading effect for the Opengl Previewer
  phong.initialize();

  // set the camera's projection matrix
  flycamera.setPerspectiveMatrix(60.0, width / (float)height, 0.1f, 100.0f);
  flycamera.setViewport(Eigen::Vector2f((float)width, (float)height));

  // load the OBJ file and materials

  Tucano::MeshImporter::loadObjFile(mesh, materials,
                                    "resources/models/dodgeColorTest.obj");


  mesh.normalizeModelMatrix();

  // pass all the materials to the Phong Shader
  for (int i = 0; i < materials.size(); ++i)
    phong.addMaterial(materials[i]);

  // set the color and size of the sphere to represent the light sources
  // same sphere is used for all sources
  lightrep.setColor(Eigen::Vector4f(1.0, 1.0, 0.0, 1.0));
  lightrep.setSize(0.15);

  // create a first ray-tracing light source at some random position
  lights.push_back(Eigen::Vector3f(-0.5, 2.0, 3.0));

  // scale the camera representation (frustum) for the ray debug
  camerarep.shapeMatrix()->scale(0.2);

  // the debug ray is a cylinder, set the radius and length of the cylinder
  ray.setSize(0.005, 1.0);

  // craete a first debug ray pointing at the center of the screen
  createDebugRay(Eigen::Vector2f(width / 2.0, height / 2.0));

  glEnable(GL_DEPTH_TEST);

  // for (int i = 0; i<mesh.getNumberOfFaces(); ++i){
  //   Tucano::Face face = mesh.getFace(i);    
  //   for (int j =0; j<face.vertex_ids.size(); ++j){
  //     std::cout<<"vid "<<j<<" "<<face.vertex_ids[j]<<std::endl;
  //     std::cout<<"vertex "<<mesh.getVertex(face.vertex_ids[j]).transpose()<<std::endl;
  //     std::cout<<"normal "<<mesh.getNormal(face.vertex_ids[j]).transpose()<<std::endl;
  //   }
  //   std::cout<<"mat id "<<face.material_id<<std::endl<<std::endl;
  //   std::cout<<"face   normal "<<face.normal.transpose() << std::endl << std::endl;
  // }



}

void Flyscene::paintGL(void) {

  // update the camera view matrix with the last mouse interactions
  flycamera.updateViewMatrix();
  Eigen::Vector4f viewport = flycamera.getViewport();

  // clear the screen and set background color
  glClearColor(0.9, 0.9, 0.9, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // position the scene light at the last ray-tracing light source
  scene_light.resetViewMatrix();
  scene_light.viewMatrix()->translate(-lights.back());

  // render the scene using OpenGL and one light source
  phong.render(mesh, flycamera, scene_light);

  // render the ray and camera representation for ray debug
  ray.render(flycamera, scene_light);
  camerarep.render(flycamera, scene_light);

  // render ray tracing light sources as yellow spheres
  for (int i = 0; i < lights.size(); ++i) {
    lightrep.resetModelMatrix();
    lightrep.modelMatrix()->translate(lights[i]);
    lightrep.render(flycamera, scene_light);
  }

  // render coordinate system at lower right corner
  flycamera.renderAtCorner();
}

void Flyscene::simulate(GLFWwindow *window) {
  // Update the camera.
  // NOTE(mickvangelderen): GLFW 3.2 has a problem on ubuntu where some key
  // events are repeated: https://github.com/glfw/glfw/issues/747. Sucks.
  float dx = (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS ? 0.05 : 0.0) -
             (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS ? 0.05 : 0.0);
  float dy = (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS
                  ? 0.05
                  : 0.0) -
             (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS
                  ? 0.05
                  : 0.0);
  float dz = (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS ? 0.05 : 0.0) -
             (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS ? 0.05 : 0.0);
  flycamera.translate(dx, dy, dz);
}

void Flyscene::createDebugRay(const Eigen::Vector2f &mouse_pos) {
  ray.resetModelMatrix();
  // from pixel position to world coordinates
  Eigen::Vector3f screen_pos = flycamera.screenToWorld(mouse_pos);

  // direction from camera center to click position
  Eigen::Vector3f dir = (screen_pos - flycamera.getCenter()).normalized();
  
	//set length
  vector<Box> boxes = getMoreBoxes();
  float closest = INFINITY;
  float hitpoint;
  Tucano::Face display_face;
  bool intersected = false;
  Eigen::Vector3f result = Eigen::Vector3f(0.9, 0.9, 0.9);
  //for each triangle check if intersects
	if (bBoxIntersection(boxes, screen_pos, flycamera.getCenter())) {
	  //std::cout << "" << std::endl;
	  for (int i = 0; i < mesh.getNumberOfFaces(); i++) {
		  //std::cout << "intersection2" << std::endl;
		  Tucano::Face current_face = mesh.getFace(i);
		  float distance;
		  if (intersect(screen_pos, flycamera.getCenter(), current_face, distance)) {
			  if (distance > 0 && closest > distance) {
				  intersected = true;
				  display_face = current_face;
				  closest = distance;

			  }
		  }
	  }
  }
	if (intersected) {
		hitpoint = closest;
		ray.setColor(Eigen::Vector4f(0, 1, 0, 1));
	}
	else {
		hitpoint = 1000;
		ray.setColor(Eigen::Vector4f(1, 0, 0, 1));
	}

  // position and orient the cylinder representing the ray
  ray.setOriginOrientation(flycamera.getCenter(), dir);
  ray.setSize(0.005, hitpoint);

  // place the camera representation (frustum) on current camera location, 
  camerarep.resetModelMatrix();
  camerarep.setModelMatrix(flycamera.getViewMatrix().inverse());
}

void Flyscene::raytraceScene(int width, int height) {
  std::cout << "ray tracing ..." << std::endl;

  // if no width or height passed, use dimensions of current viewport
  Eigen::Vector2i image_size(width, height);
  if (width == 0 || height == 0) {
    image_size = flycamera.getViewportSize();
  }

  // create 2d vector to hold pixel colors and resize to match image size
  vector<vector<Eigen::Vector3f>> pixel_data;
  pixel_data.resize(image_size[1]);
  for (int i = 0; i < image_size[1]; ++i)
    pixel_data[i].resize(image_size[0]);

  // origin of the ray is always the camera center
  Eigen::Vector3f origin = flycamera.getCenter();
  Eigen::Vector3f screen_coords;


  //make level
  int level = 0;

  vector<Box> boxes = getMoreBoxes();


  // for every pixel shoot a ray from the origin through the pixel coords
  for (int j = 0; j < image_size[1]; ++j) {
	  std::cout << j << std::endl;
    for (int i = 0; i < image_size[0]; ++i) {
      // create a ray from the camera passing through the pixel (i,j)
      screen_coords = flycamera.screenToWorld(Eigen::Vector2f(/*image_size[0] -*/ i, /*image_size[1] - */j));
      // launch raytracing for the given ray and write result to pixel data
      pixel_data[j][i] = traceRay(origin, screen_coords, boxes, level);
    }
  }

  // write the ray tracing result to a PPM image
  Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);
  std::cout << "ray tracing done! " << std::endl;
}


Eigen::Vector3f Flyscene::traceRay(Eigen::Vector3f &origin,
                                   Eigen::Vector3f &dest, vector<Box>& boxes, int level) {
	Tucano::Face display_face;
	float closest = INFINITY;
	bool intersected = false;
	Eigen::Vector3f result = Eigen::Vector3f(0.9, 0.9, 0.9);
	if (bBoxIntersection(boxes, dest, origin)) {
	for (int i = 0; i < mesh.getNumberOfFaces(); i++) {
		Tucano::Face current_face = mesh.getFace(i);
		float distance;	
			if (intersect(dest, origin, current_face, distance)) {
				if (distance > 0 && closest > distance) {
					intersected = true;
					//std::cout << "distance = " << distance << std::endl;
					display_face = current_face;
					closest = distance;
				}
			}
		}
	}
	if (intersected) result = shade(level, display_face, origin, closest * dest, boxes); // check origin*dest

	return result;
}

bool Flyscene::intersect(const Eigen::Vector3f& destination, const Eigen::Vector3f& origin, Tucano::Face& face, float& new_intersection) {

	if (face.vertex_ids.size() == 0) {
		return false;
	}
	

	Eigen::Vector3f vector_one = (mesh.getShapeModelMatrix() * mesh.getVertex(face.vertex_ids[0])).head<3>();
	Eigen::Vector3f vector_two = (mesh.getShapeModelMatrix() * mesh.getVertex(face.vertex_ids[1])).head<3>();
	Eigen::Vector3f vector_three = (mesh.getShapeModelMatrix() * mesh.getVertex(face.vertex_ids[2])).head<3>();



	// Face Normal
	Eigen::Vector3f normal = face.normal;

	float bottomTerm = (destination).dot(normal);
	if (bottomTerm == 0) {
		return false;
	}

	// Distance from the origin
	float distance = normal.dot(vector_one);

	// T-Point
	float t = (distance - origin.dot(normal)) / destination.dot(normal);

	// Intersection point
	Eigen::Vector3f intersectionPoint = origin + (t * destination);

	// Vector perpendicular to face
	Eigen::Vector3f per;

	// Edge 1
	Eigen::Vector3f edge_one = vector_two - vector_one;
	per = edge_one.cross(intersectionPoint - vector_one);
	if (normal.dot(per) < 0) {
		return false;
	}

	// Edge 2
	Eigen::Vector3f edge_two = vector_three - vector_two;
	per = edge_two.cross(intersectionPoint - vector_two);
	if (normal.dot(per) < 0) {
		return false;
	}

	// Edge 3
	Eigen::Vector3f edge_three = vector_one - vector_three;
	per = edge_three.cross(intersectionPoint - vector_three);
	if (normal.dot(per) < 0) {
		return false;
	}
	//Update the new intersection point

	new_intersection = distance3f(intersectionPoint, origin);
	return true;
}

//compute distance between 3-d vectors 
float Flyscene::distance3f(Eigen::Vector3f vec1, Eigen::Vector3f vec2) {
	return sqrt(pow(vec2(0) - vec1(0), 2) + pow(vec2(1) - vec1(1), 2) + pow(vec2(2) - vec1(2), 2));
}

// http://www.cs.utah.edu/~awilliam/box/box.pdf
bool Flyscene::bBoxIntersection(const vector<Box>& boxes, const Eigen::Vector3f& destination, const Eigen::Vector3f& origin) {

	if (destination.x() >= 0 || destination.y() >= 0 || destination.z() >= 0) {
		for (int i = 0; i < boxes.size(); i++) {
			Box box = boxes.at(i);
			float maximum_x_value, maximum_y_value, maximum_z_value, minimum_x_value, minimum_y_value, minimum_z_value;
			minimum_x_value = box.min.x();
			maximum_x_value = box.max.x();

			minimum_y_value = box.min.y();
			maximum_y_value = box.max.y();

			minimum_z_value = box.min.z();
			maximum_z_value = box.max.z();

			float minimum_tx_value, maximum_tx_value, minimum_ty_value, maximum_ty_value, minimum_tz_value, maximum_tz_value;

			// X
			if (destination.x() >= 0) {
				minimum_tx_value = (minimum_x_value - origin.x()) / destination.x();
				maximum_tx_value = (maximum_x_value - origin.x()) / destination.x();
			}
			else {
				minimum_tx_value = (maximum_x_value - origin.x()) / destination.x();
				maximum_tx_value = (minimum_x_value - origin.x()) / destination.x();
			}

			// Y
			if (destination.y() >= 0) {
				minimum_ty_value = (minimum_y_value - origin.y()) / destination.y();
				maximum_ty_value = (maximum_y_value - origin.y()) / destination.y();
			}
			else {
				minimum_ty_value = (maximum_y_value - origin.y()) / destination.y();
				maximum_ty_value = (minimum_y_value - origin.y()) / destination.y();
			}
			if ((minimum_tx_value > maximum_ty_value) || (minimum_ty_value > maximum_tx_value)) {
				continue;
			}

			if (minimum_ty_value > minimum_tx_value) {
				minimum_tx_value = minimum_ty_value;
			}
			if (maximum_ty_value < maximum_tx_value) {
				maximum_tx_value = maximum_ty_value;
			}

			// Z
			if (destination.z() >= 0) {
				minimum_tz_value = (minimum_z_value - origin.z()) / destination.z();
				maximum_tz_value = (maximum_z_value - origin.z()) / destination.z();
			}
			else {
				minimum_tz_value = (maximum_z_value - origin.z()) / destination.z();
				maximum_tz_value = (minimum_z_value - origin.z()) / destination.z();
			}

			if ((minimum_tx_value > maximum_tz_value) || (minimum_tz_value > maximum_tx_value)) {
				continue;
			}
			if (minimum_tz_value > minimum_tx_value) {
				minimum_tx_value = minimum_tz_value;
			}
			if (maximum_tz_value < maximum_tx_value) {
				maximum_tx_value = maximum_tz_value;
			}

			//std::cout << " got intersection" <<std::endl;

			return true;
		}
	}
	return false;

}

vector<Box> Flyscene::getMoreBoxes() {

	vector<Box> result;

	for (int i = 0; i < mesh.getNumberOfVertices(); i++) {

		float minX = INFINITY, minY = INFINITY, minZ = INFINITY;
		float maxX = -INFINITY, maxY = -INFINITY, maxZ = -INFINITY;
		vector<int> faceList;

		for (int j = 0; j < 1000 && (i * 1000 + j < mesh.getNumberOfVertices()); j++) {
			Eigen::Vector4f curr = /*mesh.getShapeModelMatrix() **/ mesh.getVertex(j + i * 1000);
			
			faceList.push_back(j + i*1000);

			minX = min(minX, curr.x());
			minY = min(minY, curr.y());
			minZ = min(minZ, curr.z());

			maxX = max(maxX, curr.x());
			maxY = max(maxY, curr.y());
			maxZ = max(maxZ, curr.z());
		}

		Eigen::Vector3f min = Eigen::Vector3f(minX, minY, minZ);
		Eigen::Vector3f max = Eigen::Vector3f(maxX, maxY, maxZ);

		Box resultBox = Box(min, max, faceList);
		result.push_back(resultBox);
	}

	return result;
}

Box Flyscene::getFullBox() {
	float minimum_x_value = INFINITY, minimum_y_value = INFINITY, minimum_z_value = INFINITY;
	float maximum_value_value = -INFINITY, maximum_y_value = -INFINITY, maximum_z_value = -INFINITY;
	vector<int> faceList;

	for (int i = 0; i < mesh.getNumberOfVertices(); i++) {
		Eigen::Vector4f curr = mesh.getVertex(i);
		faceList.push_back(i);

		minimum_x_value = min(minimum_x_value, curr.x());
		minimum_y_value = min(minimum_y_value, curr.y());
		minimum_z_value = min(minimum_z_value, curr.z());

		maximum_value_value = max(maximum_value_value, curr.x());
		maximum_y_value = max(maximum_y_value, curr.y());
		maximum_z_value = max(maximum_z_value, curr.z());
	}

	Eigen::Vector3f min = Eigen::Vector3f(minimum_x_value, minimum_y_value, minimum_z_value);
	Eigen::Vector3f max = Eigen::Vector3f(maximum_value_value, maximum_y_value, maximum_z_value);

	Box result = Box(min, max, faceList);
	return result;
}

Eigen::Vector3f Flyscene::computeDirectLight(Eigen::Vector4f currentColor, Tucano::Face hit, Eigen::Vector3f lightDirection, Eigen::Vector3f origin, Eigen::Vector3f destination) {

	float t;
	intersect(destination, origin, hit, t);
	Eigen::Vector3f newSmoothNormal = smoothSurfacing((origin + t * destination), hit);

	Eigen::Vector3f ka = materials[hit.material_id].getAmbient();
	Eigen::Vector3f ks = materials[hit.material_id].getSpecular();
	Eigen::Vector3f kd = materials[hit.material_id].getDiffuse();
	float shininess = materials[hit.material_id].getShininess();

	float light_intensity = currentColor.w();
	Eigen::Vector3f eyeDirection = (-1 * (destination - origin)).normalized();
	Eigen::Vector3f reflectedRay = (lightDirection - 2 * newSmoothNormal.dot(lightDirection) * newSmoothNormal).normalized();

	Eigen::Vector3f ambient = light_intensity * ka;
	Eigen::Vector3f diffuse = light_intensity * kd * std::max(lightDirection.normalized().dot(-1 * hit.normal.normalized()), 0.f);
	float cosphi = std::abs((reflectedRay).dot(eyeDirection));
	Eigen::Vector3f specular = light_intensity * ks * pow(cosphi, shininess);

	//std::cout << "ka:       " << ka.x() << " " << ka.y() << " " << ka.z() << std::endl;
	//std::cout << "ks:       " << ks.x() << " " << ks.y() << " " << ks.z() << std::endl;
	//std::cout << "kd:       " << kd.x() << " " << kd.y() << " " << kd.z() << std::endl;
	//std::cout << "ambient:  " << ambient.x() << " " << ambient.y() << " " << ambient.z() << std::endl;
	//std::cout << "diffuse:  " << diffuse.x() << " " << diffuse.y() << " " << diffuse.z() << std::endl;
	//std::cout << "cosphi:   " << cosphi << std::endl;
	//if (specular.x() > 0.1 || specular.y() > 0.1 || specular.z() > 0.1) {
	//	std::cout << "reflectedray" << reflectedRay.x() << " " << reflectedRay.y() << " " << reflectedRay.z() << std::endl;
	//	std::cout << "eyedirection" << eyeDirection.x() << " " << eyeDirection.y() << " " << eyeDirection.z() << std::endl;
	//	std::cout << "cosphi:2  " << reflectedRay.dot(eyeDirection) << std::endl;
	//	std::cout << "specular: " << specular.x() << " " << specular.y() << " " << specular.z() << std::endl;
	//}


	Eigen::Vector3f directColor = ambient + diffuse; //+specular;
	//std::cout << "color:    " << directColor.x() << " " << directColor.y() << " " << directColor.z() << std::endl;
	return directColor;
}

Eigen::Vector3f Flyscene::getCenterFace(Tucano::Face face) {
	Eigen::Vector3f v0 = mesh.getVertex(face.vertex_ids[0]).head<3>();
	Eigen::Vector3f v1 = mesh.getVertex(face.vertex_ids[1]).head<3>();
	Eigen::Vector3f v2 = mesh.getVertex(face.vertex_ids[2]).head<3>();

	Eigen::Vector3f centerFace = (v0 + v1 + v2) / 3;
	return centerFace;
}

Eigen::Vector3f Flyscene::computeReflected(Eigen::Vector3f origin, Eigen::Vector3f destination, Eigen::Vector3f lightDirection, Tucano::Face hit, int level, vector<Box>& box) {
	
	float t;
	intersect(destination, origin, hit, t);
	Eigen::Vector3f newSmoothNormal = smoothSurfacing((origin + t * destination), hit);
	
	if (level >= 1) {
		return Eigen::Vector3f(0, 0, 0);
	}
	Eigen::Vector3f normal = newSmoothNormal.head<3>();
	Eigen::Vector3f reflectedRay = (lightDirection - 2 * (lightDirection.dot(normal)) * normal).normalized(); //go through all lights

	//Ray newRay(intersection, reflectedRay);
	Eigen::Vector3f refOrigin = getCenterFace(hit);
	Eigen::Vector3f newRay = (reflectedRay - refOrigin).normalized(); //??
	Eigen::Vector3f reflectedColor = traceRay(refOrigin, newRay, box, level + 1);

	//return the total color. now reflection and refraction is 100%.
	//std::cout << "reflectedColor" << reflectedColor/* * materials[hit.material_id].getDissolveFactor()*/ << std::endl;
	return reflectedColor * materials[hit.material_id].getDissolveFactor()/* + refractedColor * materials[mat_id].getDissolveFactor()*/;
}

Eigen::Vector3f Flyscene::shade(int level, Tucano::Face hit, Eigen::Vector3f origin, Eigen::Vector3f destination, vector<Box>& box) {

	float t;
	intersect(destination, origin, hit, t);
	Eigen::Vector3f newSmoothNormal = smoothSurfacing((origin + t * destination), hit);

	if (hit.vertex_ids.size() == 0) {
		return Eigen::Vector3f(0.f, 0.f, 0.f);
	}
	Eigen::Vector4f startColor = lightrep.getColor();

	Eigen::Vector3f color = Eigen::Vector3f(0.f, 0.f, 0.f);

	for (int i = 0; i < lights.size(); i++) {
		Eigen::Vector3f lightDirection = (newSmoothNormal - lights.at(i)).normalized();
		Eigen::Vector3f directColor = computeDirectLight(startColor, hit, lightDirection, origin, destination);

		//Eigen::Vector3f reflectedRay = computeReflectedRay(hit, lightDirection, 1, 0);
		//Eigen::Vector3f refractedRay = computeRefractedRay(hit, lights.at(0).head<3>(), lightDirection);

		//float fresnelReflection = fresnelReflectedPart(lightDirection, hit);
		//float fresnelRefraction = 1 - fresnelReflection;

		//float refractedColor = 0.f;
		Eigen::Vector3f reflectedColor = computeReflected(origin, destination, lightDirection, hit, level, box);

		color = directColor + reflectedColor/* + fresnelReflection * reflectedColor + fresnelRefraction * refractedColor*/;

	}
	return color;
}

Eigen::Vector3f Flyscene::smoothSurfacing(Eigen::Vector3f hitpoint, Tucano::Face face) {

	//v0 is A
	//v1 is B
	//v2 is C

	Eigen::Vector3f v0 = (mesh.getVertex(face.vertex_ids[0]).head<3>());
	Eigen::Vector3f v1 = (mesh.getVertex(face.vertex_ids[1]).head<3>());
	Eigen::Vector3f v2 = (mesh.getVertex(face.vertex_ids[2]).head<3>());

	// Get the vertex normals
	Eigen::Vector3f normal_zero = (mesh.getNormal(face.vertex_ids[0]).head<3>());
	Eigen::Vector3f normal_one = (mesh.getNormal(face.vertex_ids[1]).head<3>());
	Eigen::Vector3f normal_two = (mesh.getNormal(face.vertex_ids[2]).head<3>());

	// Get the barycentric coordinates
	float ABCArea = ((v1 - v0).cross(v2 - v0)).size() / 2;
	float ABPArea = ((v1 - v0).cross(hitpoint - v0)).size() / 2;
	float ACPArea = ((v2 - v0).cross(hitpoint - v0)).size() / 2;
	float BCPArea = ((v1 - hitpoint).cross(v2 - hitpoint)).size() / 2;

	float u = ACPArea / ABCArea;
	float v = ABPArea / ABCArea;
	float w = BCPArea / ABCArea;

	// Calculate the new interpolated normal
	Eigen::Vector3f interpolatedNormal = w * normal_zero + u * normal_one + v * normal_two;

	// Normalize just to be sure :)
	return interpolatedNormal.normalized();
}
