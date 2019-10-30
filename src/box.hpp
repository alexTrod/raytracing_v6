#ifndef BOX_H
#define BOX_H
#include <tucano/utils/flycamera.hpp>


class Box {
public:
	Box(const Eigen::Vector3f& mi, Eigen::Vector3f& ma) {
		Box::min = mi;
		Box::max = ma;
	}

	Box createBox(Tucano::Mesh mesh) {
		for (int i = 0; i < mesh.getNumberOfVertices(); i++) {

		}
	}

	Eigen::Vector3f min;
	Eigen::Vector3f max;
};

#endif