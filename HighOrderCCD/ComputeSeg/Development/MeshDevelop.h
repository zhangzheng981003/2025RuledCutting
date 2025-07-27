#pragma once
#include "MeshLaplace.h"
#include "TriInequality.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using std::vector;
using SPMatrixXd = Eigen::SparseMatrix<double, Eigen::ColMajor>;

class MeshDevelop
{
public:
	MeshDevelop();
	void oriMesh(const Mesh& mesh);
	void tarMesh(const Mesh& mesh, double dist);
	void outMesh(Mesh& mesh);

	double develop(int maxIter = 2e4);
	
	void SetSeamStatus(const std::vector<bool>& e_status) {
		seam_status = e_status;
	}

	void GetDevelopEnergy(std::vector<double>& energy) {
		energy = dev_energy;
	}

	void SmoothStatus(bool smooth_status)
	{
		smooth_status_ = smooth_status;
	}

private:
	Mesh tmp;
	MeshTopo m_topo;
	MatrixXd m_vpfnen;
	MatrixXd m_vpfnen0;
	MatrixXd m_eb0, m_vt;
	VectorXd m_el, m_ew, m_ed, m_es, m_eA, m_fA, m_vA;
	MatrixXd m_fDhDeReU; // 9 * size

	std::vector<bool> m_fixV;
	std::vector<Eigen::JacobiSVD<MatrixXd>> m_svd;
	std::vector<Matrix3d> m_svdM, m_svdU, m_svdV;
	std::vector<Vector3d> m_svdA;

	std::vector<Eigen::SelfAdjointEigenSolver<Matrix3d>> m_eig;

	double m_wD = 1.0;
	//double m_wD = 0.1;
	double m_wP = 1e3;
	double m_similar = 0;
	double m_dist;
	double m_scale;
	double m_seam = 0.95;

	SPMatrixXd m_lap, m_elap;
	Eigen::SimplicialLDLT<SPMatrixXd> m_ldlt;
	//std::unique_ptr<AndersonAcceleration> AA;
	TriInequality triNeq;

	// seam
	std::vector<bool> seam_status;
	std::vector<double> dev_energy;
	bool smooth_status_ = false;
	double lambda_smooth_ = 0.000;
	Eigen::SimplicialLDLT<SPMatrixXd> m_s_ldlt;

private:
	void localStep();
	void globalStep();
	void updateN();
	void updateV();
	void updatePara(int iter);

	void calcEN();
	void calcFN();

	void updateD();
	void updateEU();
	void updateER();
	void updateVP();
	void updateVT();
	void updateES();

	double energy();
	double energyD();
	double energyFN();
	double energyDP();
	double energyRP();

	bool inBound();
};