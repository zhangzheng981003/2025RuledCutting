#pragma once
#include "HighOrderCCD\SplineFitting/Mesh/MeshDefinition.h"
namespace LoopGen
{
	class AngleDefectMinimizer
	{
	public:
		AngleDefectMinimizer(Mesh* m) :mesh(m) {};
		~AngleDefectMinimizer() {};
	public:
		void run();
	private:
		Mesh* mesh;
	};
}
