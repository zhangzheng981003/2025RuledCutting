#pragma once
#include <Eigen/Sparse>
#include "../SegMesh.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

class p_smooth
{
public:
	p_smooth(SegMesh& seg_mesh);

	void Run();

private:
	void SetPosition();
	void SetStepH();

	void Init();

	void SetGloToLoc();
	void SetInnerSmoothMatrix();
	void SetSeamSmoothMatrix();
	
	void Opt();

	void ComputeK();
	void UpdateEnergy();

	void ComputePartialK();
	void ComputeGradient();

	bool LinearSearch(double& alpha_k, double c2 = 0.9);


	std::vector<T> ComputeTempTriplet(int thread_id);

	template<typename Scalar>
	Scalar AngleDefect(const VH& v_h, const Eigen::Matrix<Scalar, -1, 3>& temp_pos);

private:
	SegMesh& seg_mesh_;
	Mesh& mesh_;

	Eigen::MatrixX3d ori_pos_;
	Eigen::MatrixX3d new_pos_;

	std::vector<double> step_h_;

	double step_precent_ = 0.001;

	std::vector<int> glo2inner_;
	std::vector<int> inner2glo_;
	std::vector<int> seam2glo_;

	SpMat m_sm_T_sm_;
	SpMat m_se_T_se_;

	double total_energy_;
	double smooth_energy_;
	double seam_energy_;
	double close_energy_;
	double dev_energy_;

	Eigen::VectorXd K_;
	SpMat partial_K_;

	Eigen::MatrixX3d m_gradient_;
	Eigen::MatrixX3d m_direction_;

	// Nadam
	int cout_t_ = 1;
	Eigen::MatrixX3d m_t_;
	double v_t_;
	double epsilon_e_ = 1e-8;
	double beta_1_ = 0.9;
	double beta_2_ = 0.999;

	double lambda_smooth_ = 1.0;
	double lambda_seam_ = 1.0;
	double lambda_close_ = 1.0;
	double lambda_dev_ = 1.0;

	int iter_num_ = 10000;
	int p_ = 20;

	double small_gradient_precent_ = 0.001;
	int update_lambda_num_ = 20;
	double update_lambda_precent_ = 1.1;
};

