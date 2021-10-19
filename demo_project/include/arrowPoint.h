#pragma once
#ifndef ARROWPOINT_
#define ARROWPOINT_

#include "MyMesh.h"

class ArrowPoint
{
public:
	ArrowPoint();
	ArrowPoint(MeshLib::CPoint pos) { point = pos; direction = MeshLib::CPoint(); length = 0; };
	~ArrowPoint();

	MeshLib::CPoint point;
	MeshLib::CPoint direction;
	double length;

};

ArrowPoint::ArrowPoint()
{
}

ArrowPoint::~ArrowPoint()
{
}

#endif // !ARROWPOINT_
