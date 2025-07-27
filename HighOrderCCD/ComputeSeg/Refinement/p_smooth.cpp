#include "p_smooth.h"

p_smooth::p_smooth(SegMesh& seg_mesh)
	:seg_mesh_(seg_mesh), mesh_(seg_mesh.GetMesh())
{
	SetPosition();
	SetStepH();
}

void p_smooth::Run()
{
	for (const VH& v_h : mesh_.vertices())
	{
		new_pos_.row(v_h.idx()) = Eigen::Map<Eigen::Vector3d>(mesh_.point(v_h).data(), 3);
	}

	Init();

	Opt();

	for (size_t i = 0; i < new_pos_.rows(); i++)
	{
		mesh_.set_point(VH(i),
			OpenMesh::Vec3d(new_pos_(i, 0), new_pos_(i, 1), new_pos_(i, 2)));
	}
}

void p_smooth::SetPosition()
{
	int n_v = mesh_.n_vertices();

	// position
	ori_pos_.resize(n_v, 3);
	for (const VH& v_h : mesh_.vertices())
	{
		ori_pos_.row(v_h.idx()) = Eigen::Map<Eigen::Vector3d>(mesh_.point(v_h).data(), 3);
	}

	new_pos_ = ori_pos_;
}

void p_smooth::SetStepH()
{
	step_h_.resize(mesh_.n_vertices(), 0);

	double min_len = DBL_MAX;
	double max_len = -DBL_MAX;
	for (const EH& e_h : mesh_.edges())
	{
		min_len = std::min(min_len, mesh_.calc_edge_length(e_h));
		max_len = std::max(max_len, mesh_.calc_edge_length(e_h));
	}

	for (size_t i = 0; i < mesh_.n_vertices(); i++)
	{
		step_h_[i] = step_precent_ * min_len;
	}

	lambda_dev_ *= max_len / 0.001;
}

void p_smooth::Init()
{
	std::cout << "close: " << lambda_close_
		<< " dev: " << lambda_dev_ << std::endl;

	SetGloToLoc();

	SetInnerSmoothMatrix();
	SetSeamSmoothMatrix();
}

void p_smooth::SetGloToLoc()
{
	int n_v = mesh_.n_vertices();
	const auto& seam_status = seg_mesh_.GetSeam();

	std::vector<bool> seam_v_status(n_v, false);

	HEH temp_he;
	int to_idx, from_idx;
	for (const EH& e_h : mesh_.edges())
	{
		if (!seam_status[e_h.idx()]) continue;

		temp_he = mesh_.halfedge_handle(e_h, 0);

		to_idx = mesh_.to_vertex_handle(temp_he).idx();
		from_idx = mesh_.from_vertex_handle(temp_he).idx();

		seam_v_status[to_idx] = true;
		seam_v_status[from_idx] = true;
	}

	int n_inner = 0;
	glo2inner_.assign(n_v, -1);

	seam2glo_.clear();
	inner2glo_.clear();

	for (size_t i = 0; i < n_v; i++)
	{
		if (seam_v_status[i] || mesh_.is_boundary(VH(i)))
		{
			seam2glo_.push_back(i);
		}
		else
		{
			inner2glo_.push_back(i);
			glo2inner_[i] = n_inner;
			++n_inner;
		}
	}
}

void p_smooth::SetInnerSmoothMatrix()
{
	int n_inner = inner2glo_.size();

	std::vector<T> triplet;
	triplet.reserve(n_inner * 10 * 2);
	for (size_t i = 0; i < n_inner; i++)
	{
		for (const VH& adj_v : mesh_.vv_range(VH(inner2glo_[i])))
		{
			triplet.push_back(T(i, adj_v.idx(), -1.0));
			triplet.push_back(T(i, inner2glo_[i], 1.0));
		}
	}

	SpMat m_sm(n_inner, mesh_.n_vertices());
	m_sm.setFromTriplets(triplet.begin(), triplet.end());

	m_sm_T_sm_ = m_sm.transpose() * m_sm;
}

void p_smooth::SetSeamSmoothMatrix()
{
	int n_seam = seam2glo_.size();
	const auto& seam_status = seg_mesh_.GetSeam();

	std::vector<T> triplet;
	triplet.reserve(n_seam * 10 * 2);
	for (int i = 0; i < n_seam; i++)
	{
		for (const HEH& adj_he : mesh_.voh_range(VH(seam2glo_[i])))
		{
			if (!seam_status[mesh_.edge_handle(adj_he).idx()]) continue;

			triplet.push_back(T(i, mesh_.to_vertex_handle(adj_he).idx(), -1.0));
			triplet.push_back(T(i, seam2glo_[i], 1.0));
		}
	}

	SpMat m_se(n_seam, mesh_.n_vertices());
	m_se.setFromTriplets(triplet.begin(), triplet.end());

	m_se_T_se_ = m_se.transpose() * m_se;
}

void p_smooth::Opt()
{
	const int& n_v = mesh_.n_vertices();

	double mean_len = MeshTools::AverageEdgeLength(mesh_);

	K_.resize(inner2glo_.size(), 1);
	m_gradient_.resize(new_pos_.rows(), new_pos_.cols());


	// Nadam
	v_t_ = 0;
	m_t_ = Eigen::MatrixX3d::Zero(new_pos_.rows(), new_pos_.cols());

	double prev_energy = DBL_MAX;
	double min_energy = DBL_MAX;
	Eigen::MatrixX3d min_pos = ori_pos_;

	double gamma;
	int update_cout = 0;
	for (size_t i = 0; i < iter_num_; i++)
	{
		//if (i % 1000 == 999)
		//{
		//	//lambda_dev_ *= 1.1;

		//	for (size_t i = 0; i < n_v; i++)
		//	{
		//		mesh_.set_point(VH(i),
		//			OpenMesh::Vec3d(new_pos_(i), new_pos_(i + n_v), new_pos_(i + 2 * n_v)));
		//	}

		//	MeshTools::WriteMesh(mesh_, "tmp_" + std::to_string(i) + ".obj");
		//	//MeshTools::WriteMesh(mesh_, std::to_string(i_num) + "_tmp_" + std::to_string(cont) + ".obj");

		//	std::ofstream file_out("out_" + std::to_string(i) + ".txt");
		//	//std::ofstream file_out(std::to_string(i_num) + "_out_" + std::to_string(cont) + ".txt");
		//	file_out << "cont: " << i << std::endl;
		//	file_out << "total energy: " << total_energy_
		//		<< " close energy: " << close_energy_
		//		<< " seam energy: " << seam_energy_
		//		<< " smooth energy: " << smooth_energy_
		//		<< " dev energy: " << dev_energy_ / lambda_dev_ / lambda_dev_ << std::endl;
		//}

		ComputeK();
		UpdateEnergy();

		if (min_energy > dev_energy_)
		{
			min_energy = dev_energy_;
			min_pos = new_pos_;
		}

		//if (abs(prev_energy - dev_energy_) < break_diff_)
		//{
		//	break;
		//}
		//prev_energy = dev_energy_;

		ComputePartialK();
		ComputeGradient();

		//// gradient
		//m_direction_ = -m_gradient_;

		//double max_norm = 0;
		//for (size_t i = 0; i < m_gradient_.rows(); i++)
		//{
		//	max_norm = std::max(m_gradient_.row(i).norm(), max_norm);
		//}

		//std::cout << i << " th: "
		//	<< " total energy: " << total_energy_
		//	<< " close energy: " << close_energy_
		//	<< " seam energy: " << seam_energy_
		//	<< " smooth energy: " << smooth_energy_
		//	<< " dev energy: " << dev_energy_
		//	<< " max norm: " << max_norm << std::endl;

		//bool search_status = LinearSearch(gamma, 1.0);

		//std::cout << "gamma: " << gamma
		//	<< " search status: " << search_status << std::endl;

		////system("pause");

		//if (search_status)
		//{
		//	new_pos_ = new_pos_ + gamma * m_direction_;

		//	//Projection();
		//}

		std::cout << i << " th: "
		<< " total energy: " << total_energy_
		<< " close energy: " << close_energy_
		<< " seam energy: " << seam_energy_
		<< " smooth energy: " << smooth_energy_
		<< " dev energy: " << dev_energy_ << std::endl;

		// Nadam
		m_t_ = beta_1_ * m_t_ + (1 - beta_1_) * m_gradient_;
		//m_t_ = m_t_ / (1 - pow(beta_1_, cout_t_));

		//v_t_ = beta_2_ * v_t_ + (1 - beta_2_) * (m_gradient_.transpose() * m_gradient_).trace();
		//v_t_ = v_t_ / (1 - pow(beta_2_, cout_t_));

		++cout_t_;

		m_direction_ = -m_t_;

		bool search_status = LinearSearch(gamma, 1.0);

		std::cout << "gamma: " << gamma
			<< " vt: " << v_t_
			<< " search status: " << search_status << std::endl;

		new_pos_ = new_pos_ + gamma * m_direction_;

		if (i % 100 == 99)
		{
			lambda_dev_ *= update_lambda_precent_;
		}
	}

	new_pos_ = min_pos;
}

void p_smooth::ComputeK()
{
#pragma omp parallel for
	for (int i = 0; i < inner2glo_.size(); i++)
	{
		K_(i) = AngleDefect(VH(inner2glo_[i]), new_pos_);
	}
}

void p_smooth::UpdateEnergy()
{
	close_energy_ = ((new_pos_ - ori_pos_).transpose() * (new_pos_ - ori_pos_)).trace();
	smooth_energy_ = (new_pos_.transpose() * m_sm_T_sm_ * new_pos_).trace();
	seam_energy_ = (new_pos_.transpose() * m_se_T_se_ * new_pos_).trace();
	dev_energy_ = pow(Eigen::pow(K_.array().abs(), p_).sum(), 2.0 / p_);

	total_energy_ = lambda_close_ * lambda_close_ * close_energy_
		+ lambda_smooth_ * lambda_smooth_ * smooth_energy_
		+ lambda_seam_ * lambda_seam_ * seam_energy_
		+ lambda_dev_ * lambda_dev_ * dev_energy_;

	total_energy_ *= 0.5;
}

void p_smooth::ComputePartialK()
{
	std::vector<std::vector<T>> i_triplet(6);
#pragma omp parallel for num_threads(6)
	for (int i = 0; i < 6; i++)
	{
		int thread_num = omp_get_thread_num();
		//printf("i = % d, ThreadId = % d\n", i, omp_get_thread_num());
		i_triplet[thread_num] = ComputeTempTriplet(thread_num);
	}

	std::vector<T> triplet;
	for (int i = 0; i < 6; i++)
	{
		triplet.insert(triplet.begin(), i_triplet[i].begin(), i_triplet[i].end());
	}

	partial_K_.resize(inner2glo_.size(), 3 * mesh_.n_vertices());
	partial_K_.setFromTriplets(triplet.begin(), triplet.end());
}

void p_smooth::ComputeGradient()
{
	int n_v = mesh_.n_vertices();

	m_gradient_ = lambda_close_ * lambda_close_ * (new_pos_ - ori_pos_)
		+ lambda_smooth_ * lambda_smooth_ * m_sm_T_sm_ * new_pos_
		+ lambda_seam_ * lambda_seam_ * m_se_T_se_ * new_pos_;

	//double temp_sum = pow(Eigen::pow(K_.array().abs(), p_).sum(), 2.0 / p_ - 1);

	//Eigen::VectorXd temp_vec = Eigen::pow(K_.array().abs(), p_ - 2.0).matrix().cwiseProduct(K_);

	////std::cout << "1: " << m_gradient_.rows()
	////	<< "  2: " << m_gradient_.cols()
	////	<< "  3: " << temp_vec.rows()
	////	<< "  4: " << temp_vec.cols()
	////	<< "  5: " << partial_K_.rows()
	////	<< "  6: " << partial_K_.cols()
	////	<< std::endl;

	//Eigen::VectorXd temp_gra = lambda_dev_ * lambda_dev_ * temp_sum * temp_vec.transpose() * partial_K_;
	//m_gradient_.col(0) += temp_gra.segment(0, n_v);
	//m_gradient_.col(1) += temp_gra.segment(n_v, n_v);
	//m_gradient_.col(2) += temp_gra.segment(2 * n_v, n_v);

	double p_norm = pow(Eigen::pow(K_.array().abs(), p_).sum(), 1.0 / p_);

	//std::cout << p_norm << std::endl;

	Eigen::VectorXd temp_vec_0 = K_ / p_norm;
	Eigen::VectorXd temp_vec_1 = Eigen::pow(temp_vec_0.array().abs(), p_ - 2).matrix().cwiseProduct(K_);

	Eigen::VectorXd temp_gra = lambda_dev_ * lambda_dev_ * temp_vec_1.transpose() * partial_K_;
	m_gradient_.col(0) += temp_gra.segment(0, n_v);
	m_gradient_.col(1) += temp_gra.segment(n_v, n_v);
	m_gradient_.col(2) += temp_gra.segment(2 * n_v, n_v);

	//std::cout << m_gradient_ << std::endl;
}

bool p_smooth::LinearSearch(double& alpha_i, double c2)
{
	// strong Wolfe condition
	int search_num = 200;
	double c1 = 1e-4;
	double t = 2;
	//double alpha_max = 1e10;

	double alpha_i_1 = 0;
	alpha_i = 1.0;

	double phi_0 = total_energy_;
	double phi_0_d = (m_direction_.transpose() * m_gradient_).trace();

	Eigen::MatrixX3d init_pos = new_pos_;

	bool search_status = false;
	double alpha_lo = alpha_i_1, alpha_hi = alpha_i;
	{
		double phi_i_d;

		for (size_t i = 0; i < search_num; i++)
		{
			new_pos_ = init_pos + alpha_i * m_direction_;

			//Projection();
			ComputeK();
			UpdateEnergy();

			if (total_energy_ > phi_0 + c1 * alpha_i * phi_0_d)
			{
				alpha_lo = alpha_i_1;
				alpha_hi = alpha_i;
				break;
			}

			ComputePartialK();
			ComputeGradient();

			phi_i_d = (m_direction_.transpose() * m_gradient_).trace();
			if (abs(phi_i_d) < -c2 * phi_0_d)
			{
				search_status = true;
				break;
			}

			if (phi_i_d >= 0)
			{
				alpha_lo = alpha_i;
				alpha_hi = alpha_i_1;
				break;
			}

			alpha_i_1 = alpha_i;
			alpha_i = alpha_i_1 * t;
		}
	}

	if (!search_status)
	{
		double phi_i_d;

		new_pos_ = init_pos + alpha_lo * m_direction_;
		ComputeK();
		UpdateEnergy();
		double phi_lo = total_energy_;

		int cout = 0;
		do
		{
			alpha_i = (alpha_lo + alpha_hi) / 2;

			new_pos_ = init_pos + alpha_i * m_direction_;
			ComputeK();
			UpdateEnergy();

			//std::cout << " total energy: " << total_energy_
			//	<< " close energy: " << close_energy_
			//	<< " seam energy: " << seam_energy_
			//	<< " smooth energy: " << smooth_energy_
			//	<< " dev energy: " << dev_energy_ << std::endl;

			if (total_energy_ > phi_0 + c1 * alpha_i * phi_0_d
				|| total_energy_ >= phi_lo)
			{
				alpha_hi = alpha_i;
			}
			else
			{
				ComputePartialK();
				ComputeGradient();

				phi_i_d = (m_direction_.transpose() * m_gradient_).trace();
				if (abs(phi_i_d) < -c2 * phi_0_d)
				{
					search_status = true;
					break;
				}

				if (phi_i_d * (alpha_hi - alpha_lo) >= 0)
				{
					alpha_hi = alpha_lo;
				}

				alpha_lo = alpha_i;

				new_pos_ = init_pos + alpha_lo * m_direction_;
				ComputeK();
				UpdateEnergy();
				phi_lo = total_energy_;
			}

			cout++;
		} while (cout < search_num);
	}

	new_pos_ = init_pos;
	return search_status;
}

std::vector<T> p_smooth::ComputeTempTriplet(int thread_id)
{
	int coor_idx = thread_id / 2;

	double sign;
	if (thread_id % 2 == 0)
	{
		sign = 1;
	}
	else
	{
		sign = -1;
	}

	int n_v = mesh_.n_vertices();
	std::vector<T> triplet;
	triplet.reserve(10 * n_v);

	Eigen::MatrixX3d change_pos = new_pos_;
	for (int v_id = 0; v_id < n_v; v_id++)
	{
		auto& change_row = change_pos.row(v_id);

		const double& temp_step = sign * step_h_[v_id];

		int val_id = v_id + coor_idx * n_v;

		change_row(coor_idx) += temp_step;

		if (glo2inner_[v_id] != -1)
		{
			triplet.push_back(T(
				glo2inner_[v_id],
				val_id,
				0.5 * AngleDefect(VH(v_id), change_pos) / temp_step));
		}

		for (const VH& adj_v : mesh_.vv_range(VH(v_id)))
		{
			if (glo2inner_[adj_v.idx()] != -1)
			{
				triplet.push_back(T(
					glo2inner_[adj_v.idx()],
					val_id,
					0.5 * AngleDefect(adj_v, change_pos) / temp_step));
			}
		}

		change_row(coor_idx) -= temp_step;
	}

	return triplet;
}

template<typename Scalar>
Scalar p_smooth::AngleDefect(const VH& v_h, const Eigen::Matrix<Scalar, -1, 3>& temp_pos)
{
	Scalar scalar_pi = M_PI;
	Scalar angle_defect = 2 * scalar_pi;

	int cur_idx, prev_idx, next_idx;
	Eigen::Matrix<Scalar, -1, 1> cur_vec, prev_vec;
	for (const HEH& cur_he : mesh_.voh_range(v_h))
	{
		if (!mesh_.face_handle(cur_he).is_valid()) continue;

		cur_idx = v_h.idx();
		prev_idx = mesh_.opposite_vh(cur_he).idx();
		next_idx = mesh_.to_vertex_handle(cur_he).idx();

		cur_vec = (temp_pos.row(next_idx) - temp_pos.row(cur_idx)).normalized();
		prev_vec = (temp_pos.row(cur_idx) - temp_pos.row(prev_idx)).normalized();

		angle_defect -= scalar_pi - acos(cur_vec.dot(prev_vec));
	}

	return angle_defect;
}