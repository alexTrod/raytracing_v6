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
  float dx = (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS ? 1.0 : 0.0) -
             (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS ? 1.0 : 0.0);
  float dy = (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS
                  ? 1.0
                  : 0.0) -
             (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS
                  ? 1.0
                  : 0.0);
  float dz = (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS ? 1.0 : 0.0) -
             (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS ? 1.0 : 0.0);
  flycamera.translate(dx, dy, dz);
}

void Flyscene::createDebugRay(const Eigen::Vector2f &mouse_pos) {
  ray.resetModelMatrix();
  // from pixel position to world coordinates
  Eigen::Vector3f screen_pos = flycamera.screenToWorld(mouse_pos);

  // direction from camera center to click position
  Eigen::Vector3f dir = (screen_pos - flycamera.getCenter()).normalized();
  
  // position and orient the cylinder representing the ray
  ray.setOriginOrientation(flycamera.getCenter(), dir);

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

  Box box = getFullBox();


  // for every pixel shoot a ray from the origin through the pixel coords
  for (int j = 0; j < image_size[1]; ++j) {
	  std::cout << j << std::endl;
    for (int i = 0; i < image_size[0]; ++i) {
      // create a ray from the camera passing through the pixel (i,j)
      screen_coords = flycamera.screenToWorld(Eigen::Vector2f(image_size[0] - i, image_size[1] - j));
      // launch raytracing for the given ray and write result to pixel data
      pixel_data[j][i] = traceRay(origin, screen_coords, box);
    }
  }

  // write the ray tracing result to a PPM image
  Tucano::ImageImporter::writePPMImage("result.ppm", pixel_data);
  std::cout << "ray tracing done! " << std::endl;
}


Eigen::Vector3f Flyscene::traceRay(Eigen::Vector3f &origin,
                                   Eigen::Vector3f &dest, Box& box) {
	Tucano::Face display_face;
	float closest = INFINITY;
	int level = 0;
	bool intersected = false;
	Eigen::Vector3f result = Eigen::Vector3f(0.5, 0.5, 0);
	for (int i = 0; i < mesh.getNumberOfFaces(); i++) {
		Tucano::Face current_face = mesh.getFace(i);
		float distance;
		if (bBoxIntersection(box, dest, origin)) {
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

	Eigen::Vector3f tempShading = Eigen::Vector3f(rand() / (float)RAND_MAX, rand() / (float)RAND_MAX,
		rand() / (float)RAND_MAX);
	//if (intersected) result = shade(level, display_face, origin, closest * dest); // check origin*dest
	if (intersected) result = tempShading;

	return result;
}

bool Flyscene::intersect(const Eigen::Vector3f& destination, const Eigen::Vector3f& origin, Tucano::Face& face, float& new_intersection) {

	if (face.vertex_ids.size() == 0) {
		return false;
	}

	Eigen::Vector3f vector_one = (mesh.getVertex(face.vertex_ids[0]).head<3>());
	Eigen::Vector3f vector_two = (mesh.getVertex(face.vertex_ids[1]).head<3>());
	Eigen::Vector3f vector_three = (mesh.getVertex(face.vertex_ids[2]).head<3>());



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
	//std::cout << "ffck"  << std::endl;
	//Update the new intersection point

	new_intersection = distance3f(intersectionPoint, origin);
	return true;
}

//compute distance between 3-d vectors 
float Flyscene::distance3f(Eigen::Vector3f vec1, Eigen::Vector3f vec2) {
	return sqrt(pow(vec2(0) - vec1(0), 2) + pow(vec2(1) - vec1(1), 2) + pow(vec2(2) - vec1(2), 2));
}

bool Flyscene::bBoxIntersection(const Box& box, const Eigen::Vector3f& destination, const Eigen::Vector3f& origin) {
	//std::cout << "start intersection BB" << std::endl;
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
		return false;
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
		return false;
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

Box Flyscene::getFullBox() {
	float minimum_x_value = INFINITY, minimum_y_value = INFINITY, minimum_z_value = INFINITY;
	float maximum_value_value = -INFINITY, maximum_y_value = -INFINITY, maximum_z_value = -INFINITY;

	for (int i = 0; i < mesh.getNumberOfVertices(); i++) {
		Eigen::Vector4f curr = mesh.getVertex(i);

		minimum_x_value = min(minimum_x_value, curr.x());
		minimum_y_value = min(minimum_y_value, curr.y());
		minimum_z_value = min(minimum_z_value, curr.z());

		maximum_value_value = max(maximum_value_value, curr.x());
		maximum_y_value = max(maximum_y_value, curr.y());
		maximum_z_value = max(maximum_z_value, curr.z());
	}

	Eigen::Vector3f min = Eigen::Vector3f(minimum_x_value, minimum_y_value, minimum_z_value);
	Eigen::Vector3f max = Eigen::Vector3f(maximum_value_value, maximum_y_value, maximum_z_value);

	Box result = Box(min, max);
	return result;
}
