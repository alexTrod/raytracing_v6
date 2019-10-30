#ifndef BOX_H
#define BOX_H
#include <tucano/utils/flycamera.hpp>


class Box {
public:
	Box(const Eigen::Vector3f& mi, Eigen::Vector3f& ma, vector<int> faceList) {
		Box::min = mi;
		Box::max = ma;
		Box::faceList = faceList;
	}

	int getNumberFaces() {
		return faceList.size();
	}
	int getFace(int index) {
		return faceList.at(index);
	}
	Eigen::Vector3f min;
	Eigen::Vector3f max;
	vector<int> faceList;
};

#endif