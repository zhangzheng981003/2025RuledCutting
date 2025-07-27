#include "BSplineFitting.h"
using namespace std;
BSplineFitting::BSplineFitting(Mesh& mesh) :mesh(mesh)
{
    //
}

void BSplineFitting::sample_fitting_point_on_face(int sample_num)
{
    face_point_temp.clear();
    face_point_temp.shrink_to_fit();
    face_point_temp.resize(mesh.n_faces());
    face_normal_temp.clear();
    face_normal_temp.shrink_to_fit();
    face_normal_temp.resize(mesh.n_faces());

    vector<int> ischosen(mesh.n_vertices(), 0);
    for (int i = 0; i < mesh.n_faces(); i++)
    {
        if (interested_area[i] == 1)
        {
            std::vector<Point> sample_points;
            vector<Point> v_p_3; v_p_3.clear();
            FH fh = mesh.face_handle(i);
            for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it)
            {
                Point vertex(mesh.point(*fv_it)[0], mesh.point(*fv_it)[1], mesh.point(*fv_it)[2]);
                v_p_3.push_back(vertex);

                if (ischosen[(*fv_it).idx()] == 0)
                {
                    sample_points.push_back(vertex);
                    ischosen[(*fv_it).idx()] = 1;
                }
            }
            Point A(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z()), B(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z()), C(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z());
            Triangle triangle(A, B, C);

            /*
            Point p1(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z());
            sample_points.push_back(p1);
            Point p2(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z());
            sample_points.push_back(p2);
            Point p3(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z());
            sample_points.push_back(p3);*/
            /*
            Point p1((v_p_3[0].x() + v_p_3[1].x() + v_p_3[2].x()) / 3, (v_p_3[0].y() + v_p_3[1].y() + v_p_3[2].y()) / 3, (v_p_3[0].z() + v_p_3[1].z() + v_p_3[2].z()) / 3);
            sample_points.push_back(p1);

            Point p2((v_p_3[0].x() + v_p_3[1].x() + 4 * v_p_3[2].x()) / 6, (v_p_3[0].y() + v_p_3[1].y() + 4 * v_p_3[2].y()) / 6, (v_p_3[0].z() + v_p_3[1].z() + 4 * v_p_3[2].z()) / 6);
            sample_points.push_back(p2);
            Point p3((4 * v_p_3[0].x() + v_p_3[1].x() + v_p_3[2].x()) / 6, (4 * v_p_3[0].y() + v_p_3[1].y() + v_p_3[2].y()) / 6, (4 * v_p_3[0].z() + v_p_3[1].z() + v_p_3[2].z()) / 6);
            sample_points.push_back(p3);
            Point p4((v_p_3[0].x() + 4 * v_p_3[1].x() + v_p_3[2].x()) / 6, (v_p_3[0].y() + 4 * v_p_3[1].y() + v_p_3[2].y()) / 6, (v_p_3[0].z() + 4 * v_p_3[1].z() + v_p_3[2].z()) / 6);
            sample_points.push_back(p4);

            Point p5((1 * v_p_3[0].x() + 2 * v_p_3[1].x() + 3 * v_p_3[2].x()) / 6, (1 * v_p_3[0].y() + 2 * v_p_3[1].y() + 3 * v_p_3[2].y()) / 6, (1 * v_p_3[0].z() + 2 * v_p_3[1].z() + 3 * v_p_3[2].z()) / 6);
            sample_points.push_back(p5);
            Point p6((1 * v_p_3[0].x() + 3 * v_p_3[1].x() + 2 * v_p_3[2].x()) / 6, (1 * v_p_3[0].y() + 3 * v_p_3[1].y() + 2 * v_p_3[2].y()) / 6, (1 * v_p_3[0].z() + 3 * v_p_3[1].z() + 2 * v_p_3[2].z()) / 6);
            sample_points.push_back(p6);
            Point p7((2 * v_p_3[0].x() + 1 * v_p_3[1].x() + 3 * v_p_3[2].x()) / 6, (2 * v_p_3[0].y() + 1 * v_p_3[1].y() + 3 * v_p_3[2].y()) / 6, (2 * v_p_3[0].z() + 1 * v_p_3[1].z() + 3 * v_p_3[2].z()) / 6);
            sample_points.push_back(p7);
            Point p8((2 * v_p_3[0].x() + 3 * v_p_3[1].x() + 1 * v_p_3[2].x()) / 6, (2 * v_p_3[0].y() + 3 * v_p_3[1].y() + 1 * v_p_3[2].y()) / 6, (2 * v_p_3[0].z() + 3 * v_p_3[1].z() + 1 * v_p_3[2].z()) / 6);
            sample_points.push_back(p8);
            Point p9((3 * v_p_3[0].x() + 1 * v_p_3[1].x() + 2 * v_p_3[2].x()) / 6, (3 * v_p_3[0].y() + 1 * v_p_3[1].y() + 2 * v_p_3[2].y()) / 6, (3 * v_p_3[0].z() + 1 * v_p_3[1].z() + 2 * v_p_3[2].z()) / 6);
            sample_points.push_back(p9);
            Point p10((3 * v_p_3[0].x() + 2 * v_p_3[1].x() + 1 * v_p_3[2].x()) / 6, (3 * v_p_3[0].y() + 2 * v_p_3[1].y() + 1 * v_p_3[2].y()) / 6, (3 * v_p_3[0].z() + 2 * v_p_3[1].z() + 1 * v_p_3[2].z()) / 6);
            sample_points.push_back(p10);
            Point p11((1 * v_p_3[0].x() + 1 * v_p_3[1].x() + 3 * v_p_3[2].x()) / 5, (1 * v_p_3[0].y() + 1 * v_p_3[1].y() + 3 * v_p_3[2].y()) / 5, (1 * v_p_3[0].z() + 1 * v_p_3[1].z() + 3 * v_p_3[2].z()) / 5);
            sample_points.push_back(p11);
            Point p12((1 * v_p_3[0].x() + 3 * v_p_3[1].x() + 1 * v_p_3[2].x()) / 5, (1 * v_p_3[0].y() + 3 * v_p_3[1].y() + 1 * v_p_3[2].y()) / 5, (1 * v_p_3[0].z() + 3 * v_p_3[1].z() + 1 * v_p_3[2].z()) / 5);
            sample_points.push_back(p12);
            Point p13((3 * v_p_3[0].x() + 1 * v_p_3[1].x() + 1 * v_p_3[2].x()) / 5, (3 * v_p_3[0].y() + 1 * v_p_3[1].y() + 1 * v_p_3[2].y()) / 5, (3 * v_p_3[0].z() + 1 * v_p_3[1].z() + 1 * v_p_3[2].z()) / 5);
            sample_points.push_back(p13);
            Point p14((1 * v_p_3[0].x() + 2 * v_p_3[1].x() + 2 * v_p_3[2].x()) / 5, (1 * v_p_3[0].y() + 2 * v_p_3[1].y() + 2 * v_p_3[2].y()) / 5, (1 * v_p_3[0].z() + 2 * v_p_3[1].z() + 2 * v_p_3[2].z()) / 5);
            sample_points.push_back(p14);
            Point p15((2 * v_p_3[0].x() + 1 * v_p_3[1].x() + 2 * v_p_3[2].x()) / 5, (2 * v_p_3[0].y() + 1 * v_p_3[1].y() + 2 * v_p_3[2].y()) / 5, (2 * v_p_3[0].z() + 1 * v_p_3[1].z() + 2 * v_p_3[2].z()) / 5);
            sample_points.push_back(p15);
            Point p16((2 * v_p_3[0].x() + 2 * v_p_3[1].x() + 1 * v_p_3[2].x()) / 5, (2 * v_p_3[0].y() + 2 * v_p_3[1].y() + 1 * v_p_3[2].y()) / 5, (2 * v_p_3[0].z() + 2 * v_p_3[1].z() + 1 * v_p_3[2].z()) / 5);
            sample_points.push_back(p16);
            Point p17((1 * v_p_3[0].x() + 1 * v_p_3[1].x() + 2 * v_p_3[2].x()) / 4, (1 * v_p_3[0].y() + 1 * v_p_3[1].y() + 2 * v_p_3[2].y()) / 4, (1 * v_p_3[0].z() + 1 * v_p_3[1].z() + 2 * v_p_3[2].z()) / 4);
            sample_points.push_back(p17);
            Point p18((1 * v_p_3[0].x() + 2 * v_p_3[1].x() + 1 * v_p_3[2].x()) / 4, (1 * v_p_3[0].y() + 2 * v_p_3[1].y() + 1 * v_p_3[2].y()) / 4, (1 * v_p_3[0].z() + 2 * v_p_3[1].z() + 1 * v_p_3[2].z()) / 4);
            sample_points.push_back(p18);
            Point p19((2 * v_p_3[0].x() + 1 * v_p_3[1].x() + 1 * v_p_3[2].x()) / 4, (2 * v_p_3[0].y() + 1 * v_p_3[1].y() + 1 * v_p_3[2].y()) / 4, (2 * v_p_3[0].z() + 1 * v_p_3[1].z() + 1 * v_p_3[2].z()) / 4);
            sample_points.push_back(p19);
            */
            std::vector<Vec3> sample_points_vec3; sample_points_vec3.clear();
            for (unsigned int j = 0; j < sample_points.size(); ++j)
            {
                Vec3 p;
                p[0] = sample_points[j].x();
                p[1] = sample_points[j].y();
                p[2] = sample_points[j].z();
                sample_points_vec3.push_back(p);
            }
            face_point_temp[i] = sample_points_vec3;
        }
    }
}

void BSplineFitting::BSpline_initialization()
{
    //1.p2p
    sample_fitting_point_on_face(20);
    cout << mesh2surface_temp.size() << endl;
    //cover_area_temp.clear();
    int siz = mesh2surface_temp.size();
    for (int i = 0; i < siz; i++)
    {
        mesh2surface_temp.pop_back();
    }
    //vector<Mesh2Surface>().swap(mesh2surface_temp);
    cout << mesh2surface_temp.size() << endl;
    mesh2surface_temp.shrink_to_fit();
    for (int i = 0; i < cover_area_temp.size(); i++)
    {
        for (int j = 0; j < face_point_temp[cover_area_temp[i]].size(); j++)
        {
            Mesh2Surface m2f;
            m2f.point_on_mesh[0] = face_point_temp[cover_area_temp[i]][j].x();
            m2f.point_on_mesh[1] = face_point_temp[cover_area_temp[i]][j].y();
            m2f.point_on_mesh[2] = face_point_temp[cover_area_temp[i]][j].z();

            m2f.normal = is_neg * mesh.calc_face_normal(mesh.face_handle(cover_area_temp[i]));


            mesh2surface_temp.push_back(m2f);
        }
    }

    vector<Standard_Real> u_sample = generate_equally_spaced_vector(100, interval_up_bound);
    vector<Standard_Real> v_sample = generate_equally_spaced_vector(20, 1);
    my_cgal_triangles= std::vector<MyTriangle>();
    my_cgal_triangles.shrink_to_fit();
    cgal_triangles= std::vector<Triangle>();
    cgal_triangles.shrink_to_fit();
    for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
    {
        vector<Point> v_p_3;
        v_p_3.clear();
        for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
        {
            Point vertex(mesh.point(*fv_it)[0], mesh.point(*fv_it)[1], mesh.point(*fv_it)[2]);
            v_p_3.push_back(vertex);
        }

        auto f_normal = mesh.calc_face_normal(*f_it);
        Vec3 f_n;
        f_n[0] = f_normal[0];
        f_n[1] = f_normal[1];
        f_n[2] = f_normal[2];
        my_cgal_triangles.emplace_back(
            Point(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z()),
            Point(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z()),
            Point(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z()), f_it->idx(), f_n);
        cgal_triangles.emplace_back(
            Point(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z()),
            Point(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z()),
            Point(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z()));
        v_p_3.clear();
    }

    my_tree.clear();
    tree.clear();
    my_tree.insert(my_cgal_triangles.begin(), my_cgal_triangles.end());
    my_tree.build();
    my_tree.accelerate_distance_queries();

    tree.insert(cgal_triangles.begin(), cgal_triangles.end());
    tree.build();
    tree.accelerate_distance_queries();

}

std::vector<Standard_Real> BSplineFitting::generate_equally_spaced_vector(int num_elements, double up_bound)
{
    if (num_elements <= 1)
    {
        std::cerr << "Invalid input. Number of elements should be a positive integer." << std::endl;
        return std::vector<Standard_Real>();
    }

    std::vector<Standard_Real> result;
    Standard_Real step = up_bound / (num_elements - 1); // 计算等分的步长

    for (int i = 0; i < num_elements; ++i)
    {
        Standard_Real value = i * step; // 计算等分的值
        result.push_back(value);
    }

    return result;
}

void BSplineFitting::create_BSpline_curve(vector<CtrInfo>& ctr_point, Handle(Geom_BSplineCurve)& curve)
{
    int ctr_num = ctr_point.size();
    TColgp_Array1OfPnt Poles(1, ctr_num);
    TColStd_Array1OfReal PolesWeight(1, ctr_num);
    for (int i = 0; i < ctr_num; i++)
    {
        Poles.SetValue(i + 1, ctr_point[i].position);
        PolesWeight.SetValue(i + 1, ctr_point[i].weight);
    }

    Standard_Integer PNum = ctr_point.size();

    Standard_Integer KNum = PNum - poly_degree + 1;

    TColStd_Array1OfReal knots(1, KNum);
    TColStd_Array1OfInteger mults(1, KNum);
    for (int i = 0; i < KNum; ++i)
    {
        //cout << i <<" " << u_knot_temp[i + poly_degree] << "success" << endl;
        knots.SetValue(i + 1, u_knot_temp[i + poly_degree]);

        if (i == 0 || i == KNum - 1)
        {
            mults.SetValue(i + 1, poly_degree + 1);
        }
        else
        {
            mults.SetValue(i + 1, 1);
        }

    }

    curve = new Geom_BSplineCurve(Poles, PolesWeight, knots, mults, poly_degree);

}

void BSplineFitting::create_BSpline_surface(vector<CtrInfo>& ctr_point1, vector<CtrInfo>& ctr_point2, Handle(Geom_BSplineSurface)& surface)
{
    Handle(Geom_BSplineCurve) curve1 = nullptr;
    create_BSpline_curve(ctr_point1, curve1);
    Handle(Geom_BSplineCurve) curve2 = nullptr;
    create_BSpline_curve(ctr_point2, curve2);

    GeomFill_FillingStyle Type = GeomFill_StretchStyle;
    GeomFill_BSplineCurves aGeomFill(curve1, curve2, Type);
    surface = aGeomFill.Surface();
}

void BSplineFitting::BSpline_surface_viewer(Handle(Geom_BSplineSurface)& surface, int sample_num)
{
    vector<Standard_Real> u_temp = generate_equally_spaced_vector(sample_num, interval_up_bound);
    vector<Standard_Real> v_temp = generate_equally_spaced_vector(2, 1.0);

    //2.Bezier_info构造
    Matrix<BSplineFitting::BSplinePoint, Eigen::Dynamic, Eigen::Dynamic> bsp_temp;
    bsp_temp.resize(u_temp.size(), v_temp.size());
    for (int i = 0; i < u_temp.size(); i++)
    {
        for (int j = 0; j < v_temp.size(); j++)
        {
            BSplinePoint bsp;
            bsp.u = u_temp[i];
            bsp.v = v_temp[j];
            surface->D0(bsp.u, bsp.v, bsp.position);
            bsp.bsp_position[0] = bsp.position.X();
            bsp.bsp_position[1] = bsp.position.Y();
            bsp.bsp_position[2] = bsp.position.Z();
            bsp_temp(i, j) = bsp;
        }
    }

    //3.
    Mesh outmesh;
    Matrix<int, Eigen::Dynamic, Eigen::Dynamic> v_idx_mat;
    v_idx_mat.resize(u_temp.size(), v_temp.size());
    int init_v_idx = 0;
    for (int i = 0; i < bsp_temp.rows(); i++)
    {
        for (int j = 0; j < bsp_temp.cols(); j++)
        {
            Mesh::VertexHandle v0 = outmesh.add_vertex(Mesh::Point(bsp_temp(i, j).bsp_position.data()));
            v_idx_mat(i, j) = init_v_idx;
            init_v_idx++;
        }
    }

    for (int i = 0; i < u_temp.size() - 1; ++i)
    {
        for (int j = 0; j < v_temp.size() - 1; ++j)
        {
            Mesh::VertexHandle v0 = outmesh.vertex_handle(v_idx_mat(i, j));
            Mesh::VertexHandle v1 = outmesh.vertex_handle(v_idx_mat(i + 1, j));
            Mesh::VertexHandle v2 = outmesh.vertex_handle(v_idx_mat(i + 1, j + 1));
            Mesh::VertexHandle v3 = outmesh.vertex_handle(v_idx_mat(i, j + 1));

            std::vector<Mesh::VertexHandle> face_vhandles0;
            face_vhandles0.push_back(v0);
            face_vhandles0.push_back(v1);
            face_vhandles0.push_back(v2);
            outmesh.add_face(face_vhandles0);

            std::vector<Mesh::VertexHandle> face_vhandles1;
            face_vhandles1.push_back(v0);
            face_vhandles1.push_back(v2);
            face_vhandles1.push_back(v3);
            outmesh.add_face(face_vhandles1);
        }
    }

    //4.构造并输出网格
    MeshTools::WriteMesh(mesh, "./inputmesh.obj");
    MeshTools::WriteMesh(outmesh, "./outputmesh_" + to_string(0) + "_" + to_string(0) + ".obj");
}

void BSplineFitting::BSpline_surface_viewer_2(const Data& spline, int u_sample_num, int v_sample_num, int out_number, int it_num)
{
    VectorXd BSpline;
    BSpline.resize(spline.rows() * 3);
    for (int i = 0; i < spline.rows(); i++)
    {

        BSpline[3 * i] = spline.coeff(i, 0);
        BSpline[3 * i + 1] = spline.coeff(i, 1);
        BSpline[3 * i + 2] = spline.coeff(i, 2);
    }


    vector<vector<CtrInfo>> ctr_point_temp;
    int curve_ctr_num = BSpline.size() / 6;

    //for (int i = 0; i < 2; i++)
    //{
    vector<CtrInfo> ctr_point1;
    for (int j = 0; j < curve_ctr_num; j++)
    {
        CtrInfo ci;
        Standard_Real px = BSpline[3 * j];
        Standard_Real py = BSpline[3 * j + 1];
        Standard_Real pz = BSpline[3 * j + 2];
        ci.position = { px,py,pz };
        ci.weight = 1;
        ctr_point1.push_back(ci);
    }
    ctr_point_temp.push_back(ctr_point1);

    vector<CtrInfo> ctr_point2;
    for (int j = 0; j < curve_ctr_num; j++)
    {
        CtrInfo ci;
        Standard_Real px = BSpline[3 * (j + curve_ctr_num)];
        Standard_Real py = BSpline[3 * (j + curve_ctr_num) + 1];
        Standard_Real pz = BSpline[3 * (j + curve_ctr_num) + 2];
        ci.position = { px,py,pz };
        ci.weight = 1;
        ctr_point2.push_back(ci);
    }
    ctr_point_temp.push_back(ctr_point2);
    //}

    Handle(Geom_BSplineSurface) surface = nullptr;
    create_BSpline_surface(ctr_point_temp[0], ctr_point_temp[1], surface);
    //interval_up_bound = divide_u_knot[2][0].back();
    interval_up_bound = u_knot_temp.back();
    //cout << interval_up_bound << endl;
    vector<Standard_Real> u_temp;// = generate_equally_spaced_vector(sample_num, interval_up_bound);

    u_temp.resize((u_knot_temp.size() - 3 * 2 - 1) * u_sample_num + 1);
    for (int i = 0; i < (u_knot_temp.size() - poly_degree * 2 - 1); i++)
    {
        for (int j = 0; j < u_sample_num; j++)
        {
            u_temp[i * u_sample_num + j] = (u_knot_temp[i + 1 + poly_degree] * (j)+u_knot_temp[i + poly_degree] * (u_sample_num - j)) / u_sample_num;
        }
    }
    u_temp[(u_knot_temp.size() - 3 * 2 - 1) * u_sample_num] = u_knot_temp.back();
    vector<Standard_Real> v_temp = generate_equally_spaced_vector(v_sample_num, 1);

    //2.Bezier_info构造
    Matrix<BSplineFitting::BSplinePoint, Eigen::Dynamic, Eigen::Dynamic> bsp_temp;
    bsp_temp.resize(u_temp.size(), v_temp.size());
    for (int i = 0; i < u_temp.size(); i++)
    {
        for (int j = 0; j < v_temp.size(); j++)
        {
            BSplinePoint bsp;
            bsp.u = u_temp[i];
            bsp.v = v_temp[j];
            surface->D0(bsp.u, bsp.v, bsp.position);
            bsp.bsp_position[0] = bsp.position.X();
            bsp.bsp_position[1] = bsp.position.Y();
            bsp.bsp_position[2] = bsp.position.Z();
            bsp_temp(i, j) = bsp;
        }
    }

    //3.
    Mesh outmesh;
    Matrix<int, Eigen::Dynamic, Eigen::Dynamic> v_idx_mat;
    v_idx_mat.resize(u_temp.size(), v_temp.size());
    int init_v_idx = 0;
    for (int i = 0; i < bsp_temp.rows(); i++)
    {
        for (int j = 0; j < bsp_temp.cols(); j++)
        {
            Mesh::VertexHandle v0 = outmesh.add_vertex(Mesh::Point(bsp_temp(i, j).bsp_position.data()));
            v_idx_mat(i, j) = init_v_idx;
            init_v_idx++;
        }
    }

    for (int i = 0; i < u_temp.size() - 1; ++i)
    {
        for (int j = 0; j < v_temp.size() - 1; ++j)
        {
            Mesh::VertexHandle v0 = outmesh.vertex_handle(v_idx_mat(i, j));
            Mesh::VertexHandle v1 = outmesh.vertex_handle(v_idx_mat(i + 1, j));
            Mesh::VertexHandle v2 = outmesh.vertex_handle(v_idx_mat(i + 1, j + 1));
            Mesh::VertexHandle v3 = outmesh.vertex_handle(v_idx_mat(i, j + 1));

            std::vector<Mesh::VertexHandle> face_vhandles0;
            face_vhandles0.push_back(v0);
            face_vhandles0.push_back(v1);
            face_vhandles0.push_back(v2);
            outmesh.add_face(face_vhandles0);

            std::vector<Mesh::VertexHandle> face_vhandles1;
            face_vhandles1.push_back(v0);
            face_vhandles1.push_back(v2);
            face_vhandles1.push_back(v3);
            outmesh.add_face(face_vhandles1);
        }
    }

    //4.构造并输出网格
    MeshTools::WriteMesh(mesh, "./inputmesh.obj");
    if (out_number == -2)
    {
        MeshTools::WriteMesh(outmesh, "D:/programfiles/3/11/build/result/" + result_path_fitting1 + "/123/" + to_string(it_num) + "/res.obj");
    }
    else
    {
        MeshTools::WriteMesh(outmesh, "D:/programfiles/3/11/build/result/" + result_path_fitting1 + "/123/" + to_string(it_num) + "/outputmesh_" + to_string(out_number) + ".obj");
    }
    std::cout << "success" << endl << endl;
}

void BSplineFitting::run()
{
    vector<vector<CtrInfo>> ctr_point_temp;
    int ctr_num = 9;
    for (int i = 0; i < 2; i++)
    {
        vector<CtrInfo> ctr_point;
        for (int k = 0; k < ctr_num; k++)
        {
            CtrInfo ci;
            ci.position = { Standard_Real(i),Standard_Real(k),Standard_Real(0) };
            ci.weight = 1;
            ctr_point.push_back(ci);
        }
        ctr_point_temp.push_back(ctr_point);
    }

    Handle(Geom_BSplineSurface) surface = nullptr;
    create_BSpline_surface(ctr_point_temp[0], ctr_point_temp[1], surface);
    BSpline_surface_viewer(surface, 30);
}

int BSplineFitting::u_k_position(double u, vector<double>& u_knot)
{
    for (int i = 0; i < u_knot.size(); i++)
    {
        if (u == u_knot[0])
        {
            return 0;
        }
        else if (u < u_knot[i])
        {
            return i - 1;
        }
        else if (u == u_knot.back())
        {
            return u_knot.size() - 1;
        }
    }
}

void BSplineFitting::calc_basis_fuction(double u, int k, int poly_degree, std::vector<double>& basis_func)
{
    /*!
    *\brief b样条曲线基函数计算
    *\ param double u 参数值
    *\ param int k 区间 k。可以通过std::upper_bound()-1获得
    *\ param std::vector<double> & basis_func basis_func基函数值，对应N_{k-p},...,N_k
    *\ Returns:   void
    */

    const int& p = poly_degree;
    const std::vector<double>& knots = u_knot_temp;
    basis_func.clear();
    basis_func.shrink_to_fit();
    basis_func.resize(p + 1);

    //2p+1个 N_{i,0}
    int n = 2 * p + 1;
    vector<double> temp(n, 0);
    temp[p] = 1;

    //迭代p次
    for (int j = 1; j <= p; ++j)
    {
        //区间 [k-p,k+p+1)
        for (int i = k - p, h = 0; h < (n - j); ++h, ++i)
        {
            //递推公式
            double a = (u - knots[i]);
            double dev = (knots[i + j] - knots[i]);
            a = (dev != 0) ? a / dev : 0;

            double b = (knots[i + j + 1] - u);
            dev = (knots[i + j + 1] - knots[i + 1]);
            b = (dev != 0) ? b / dev : 0;

            temp[h] = a * temp[h] + b * temp[h + 1];
        }
    }

    //拷贝前 p+1个值到basis_func
    std::copy(temp.begin(), temp.begin() + p + 1, basis_func.begin());
}

// 优化
void BSplineFitting::calc_data_term_grad_hessien(Eigen::VectorXd& BSpline, scalar_t& data_value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian)
{
    ofstream ioszz1("touying" + to_string(iter_num) + ".txt");
    data_value = 0;
    int ctr_dim = BSpline.size();//X维数
    grad.resize(ctr_dim);
    grad.setZero();
    hessian.resize(ctr_dim, ctr_dim);
    hessian.setZero();

    vector<VectorXd> m2s_grad_temp; m2s_grad_temp.clear();
    vector<MatrixXd> m2s_hessian_temp; m2s_hessian_temp.clear();
    vector<double> m2s_dist; m2s_dist.clear();

    m2s_weight.clear();
    m2s_weight.shrink_to_fit();

    //初始化自变量
    //1.
    Vec12 X;
    X.resize(ctr_dim);
    for (int r = 0; r < ctr_dim / 3; r++)
    {
        X(3 * r).value() = BSpline(3 * r);
        X(3 * r + 1).value() = BSpline(3 * r + 1);
        X(3 * r + 2).value() = BSpline(3 * r + 2);
    }

    // 2.初始化一阶导数; repeat partial derivatives for the inner AutoDiffScalar
    for (int id = 0; id < ctr_dim; id++)
    {
        X(id).derivatives().resize(ctr_dim);
        X(id).derivatives().setZero();
        X(id).derivatives()(id) = 1;
        X(id).value().derivatives() = inner_derivative_t::Unit(ctr_dim, id);
    }

    //  3.初始化海塞矩阵;set the hessian matrix to zero
    for (int idx = 0; idx < ctr_dim; idx++)
    {
        for (int id = 0; id < ctr_dim; id++)
        {
            X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(ctr_dim);
        }
    }

    std::vector<vector<Vec3>> P;
    int curve_ctr_num = ctr_dim / 6;
    P.resize(curve_ctr_num);
    for (int i = 0; i < P.size(); i++)
    {
        P[i].resize(2);
    }

    for (int i = 0; i < P.size(); i++)
    {
        P[i][0][0] = X[3 * i];
        P[i][0][1] = X[3 * i + 1];
        P[i][0][2] = X[3 * i + 2];
        P[i][1][0] = X[3 * curve_ctr_num + 3 * i];
        P[i][1][1] = X[3 * curve_ctr_num + 3 * i + 1];
        P[i][1][2] = X[3 * curve_ctr_num + 3 * i + 2];
    }

    //mesh2surface
    scalar_t data_mesh2surface = 0;
    MatrixXd data_mesh2surface_hessian;
    VectorXd data_mesh2surface_grad;
    data_mesh2surface_grad.resize(ctr_dim);
    data_mesh2surface_grad.setZero();
    data_mesh2surface_hessian.resize(ctr_dim, ctr_dim);
    data_mesh2surface_hessian.setZero();

    //int curve_ctr_num = BSpline.size() / 6;
    vector<CtrInfo> ctr_point1;
    for (int j = 0; j < curve_ctr_num; j++)
    {
        CtrInfo ci;
        Standard_Real px = BSpline[3 * j];
        Standard_Real py = BSpline[3 * j + 1];
        Standard_Real pz = BSpline[3 * j + 2];
        ci.position = { px,py,pz };
        ci.weight = 1;
        ctr_point1.push_back(ci);
    }

    vector<CtrInfo> ctr_point2;
    for (int j = 0; j < curve_ctr_num; j++)
    {
        CtrInfo ci;
        Standard_Real px = BSpline[3 * (j + curve_ctr_num)];
        Standard_Real py = BSpline[3 * (j + curve_ctr_num) + 1];
        Standard_Real pz = BSpline[3 * (j + curve_ctr_num) + 2];
        ci.position = { px,py,pz };
        ci.weight = 1;
        ctr_point2.push_back(ci);
    }

    Handle(Geom_BSplineSurface) BSpline_surface = nullptr;
    create_BSpline_surface(ctr_point1, ctr_point2, BSpline_surface);
    double nor_sum = 0;
    int false_num = 0;
    for (int k = 0; k < mesh2surface_temp.size(); k++)
    {
        if (1)
        {
            gp_Pnt pointToProject(mesh2surface_temp[k].point_on_mesh[0].value().value(), mesh2surface_temp[k].point_on_mesh[1].value().value(), mesh2surface_temp[k].point_on_mesh[2].value().value());
            GeomAPI_ProjectPointOnSurf projectPoint(pointToProject, BSpline_surface);
            
            if (projectPoint.NbPoints() > 0)
            {
                double dot_temp = 0.0;
                int idx_max = 0;



                int flag = 0;
                double e_temp = 0;
                VectorXd grad_temp;

                Eigen::MatrixXd B_temp;
                B_temp.resize(ctr_dim, ctr_dim);
                vector<double> basis_temp;
                Standard_Real u_temp, v_temp;
                int u_interval_order_temp;
                Standard_Real dist_min = 1e6;
                Vec3 bspline_surface_value1;
                bspline_surface_value1.setZero();
                //cout << "投影点个数为： " << projectPoint.NbPoints() << endl << endl;
                scalar_t e_temp1;
                for (int idx = 1; idx <= projectPoint.NbPoints(); idx++)
                {
                    mesh2surface_temp[k].is_found = true;
                    Standard_Real u, v;
                    projectPoint.Parameters(idx, u, v);



                    Standard_Real dist;
                    dist = projectPoint.Distance(idx);




                    // if(projectPoint.NbPoints()>1) cout << "投影点个数" << projectPoint.NbPoints() << endl;
                flag1:mesh2surface_temp[k].u = u;
                    mesh2surface_temp[k].v = v;

                    //cout << u << " " << v << endl;

                    //auto it_u = std::upper_bound(u_knot_temp.begin(), u_knot_temp.end(), mesh2surface_temp[k].u);
                    int u_interval_order = u_k_position(mesh2surface_temp[k].u, u_knot_temp);


                    vector<double> u_basis_temp;
                    if (u_interval_order == 0)
                    {
                        vector<double> basis_temp(poly_degree + 1, 0);
                        basis_temp[0] = 1;
                        u_basis_temp = basis_temp;
                    }
                    else if (u_interval_order == u_knot_temp.size() - 1)
                    {
                        vector<double> basis_temp(poly_degree + 1, 0);
                        basis_temp[poly_degree] = 1;
                        u_basis_temp = basis_temp;
                    }
                    else
                    {
                        calc_basis_fuction(mesh2surface_temp[k].u, u_interval_order, poly_degree, u_basis_temp);
                    }

                    // cout << mesh2surface_temp[k].u << "  " << u_interval_order << endl;
                    mesh2surface_temp[k].u_basis_temp = u_basis_temp;
                    mesh2surface_temp[k].u_interval_order = u_interval_order;

                    Vec3 bspline_surface_value;
                    bspline_surface_value.setZero();

                    if (u_interval_order == 0)
                    {
                        Vec3 pi0, pi1;
                        pi0 = P[0][0];
                        pi1 = P[0][1];
                        bspline_surface_value[0] = mesh2surface_temp[k].v * pi1[0] + (1 - mesh2surface_temp[k].v) * pi0[0];
                        bspline_surface_value[1] = mesh2surface_temp[k].v * pi1[1] + (1 - mesh2surface_temp[k].v) * pi0[1];
                        bspline_surface_value[2] = mesh2surface_temp[k].v * pi1[2] + (1 - mesh2surface_temp[k].v) * pi0[2];
                    }
                    else if (u_interval_order == u_knot_temp.size() - 1)
                    {
                        Vec3 pi0, pi1;
                        pi0 = P[curve_ctr_num - 1][0];
                        pi1 = P[curve_ctr_num - 1][1];
                        bspline_surface_value[0] = mesh2surface_temp[k].v * pi1[0] + (1 - mesh2surface_temp[k].v) * pi0[0];
                        bspline_surface_value[1] = mesh2surface_temp[k].v * pi1[1] + (1 - mesh2surface_temp[k].v) * pi0[1];
                        bspline_surface_value[2] = mesh2surface_temp[k].v * pi1[2] + (1 - mesh2surface_temp[k].v) * pi0[2];
                    }
                    else
                    {
                        for (int i = 0; i < u_basis_temp.size(); i++)
                        {
                            Vec3 pi0, pi1;
                            pi0 = P[u_interval_order - poly_degree + i][0];
                            pi1 = P[u_interval_order - poly_degree + i][1];
                            bspline_surface_value[0] += u_basis_temp[i] * mesh2surface_temp[k].v * pi1[0] + u_basis_temp[i] * (1 - mesh2surface_temp[k].v) * pi0[0];
                            bspline_surface_value[1] += u_basis_temp[i] * mesh2surface_temp[k].v * pi1[1] + u_basis_temp[i] * (1 - mesh2surface_temp[k].v) * pi0[1];
                            bspline_surface_value[2] += u_basis_temp[i] * mesh2surface_temp[k].v * pi1[2] + u_basis_temp[i] * (1 - mesh2surface_temp[k].v) * pi0[2];
                        }
                    }

                    scalar_t e = (-bspline_surface_value + mesh2surface_temp[k].point_on_mesh) * (-bspline_surface_value + mesh2surface_temp[k].point_on_mesh).transpose();

                    //data_mesh2surface += e;

                    double dot_nor = 0;
                    auto temp1 = mesh2surface_temp[k].normal;
                    temp1[0] = (bspline_surface_value - mesh2surface_temp[k].point_on_mesh)[0].value().value();
                    temp1[1] = (bspline_surface_value - mesh2surface_temp[k].point_on_mesh)[1].value().value();
                    temp1[2] = (bspline_surface_value - mesh2surface_temp[k].point_on_mesh)[2].value().value();
                    temp1.normalize();
                    if (iter_num == 0)
                    {
                        dot_nor = temp1 | mesh2surface_temp[k].normal;
                        nor_sum += dot_nor;
                    }
                    else
                    {
                        dot_nor = temp1 | mesh2surface_temp[k].normal;
                    }




                    //data_mesh2surface_grad += e.value().derivatives();
                    dist_min = projectPoint.LowerDistance();
                    Eigen::MatrixXd B;
                    B.resize(ctr_dim, ctr_dim);

                    for (int r = 0; r < ctr_dim; r++)
                    {
                        B.row(r) = e.derivatives()(r).derivatives().transpose();
                    }
                    //cout << "u= " << mesh2surface_temp[k].u << " " << u_interval_order << endl;
                    //data_mesh2surface_hessian += B;

                    //cout <<"指标： "<<idx <<endl<< "内积大小: " << dot_nor << endl << "距离为： " << dist << endl<<endl<<endl;
                    if (projectPoint.NbPoints() == 1 && (iter_num == 0 || dot_nor > 0.0))
                    {
                        //dist_min = dist;
                        bspline_surface_value1 = bspline_surface_value;

                        dot_temp = dot_nor;
                        idx_max = idx;
                        B_temp = B;
                        e_temp = e.value().value();
                        e_temp1 = e;
                        grad_temp = e.value().derivatives();
                        basis_temp = u_basis_temp;
                        u_temp = u;
                        v_temp = v;


                        u_interval_order_temp = u_interval_order;
                        idx = projectPoint.NbPoints();
                        flag = 1;
                    }
                    else if (dist < 1.3 * dist_min && ((dot_nor > 1.00 * dot_temp) || (iter_num == 0 && abs(dot_nor) > abs(dot_temp))) && u >= u_knot_temp.front() && u <= u_knot_temp.back())
                    {


                        e_temp1 = e;
                        bspline_surface_value1 = bspline_surface_value;
                        //ioszz1 << mesh2surface_temp[k].point_on_mesh[0] << " " << mesh2surface_temp[k].point_on_mesh[1] << " " << mesh2surface_temp[k].point_on_mesh[2] << " " << bspline_surface_value1[0].value().value() << " " << bspline_surface_value1[1].value().value() << " " << bspline_surface_value1[2].value().value() << endl;
                        dot_temp = dot_nor;
                        idx_max = idx;
                        B_temp = B;
                        e_temp = e.value().value();
                        grad_temp = e.value().derivatives();
                        basis_temp = u_basis_temp;
                        u_temp = u;
                        v_temp = v;
                        u_interval_order_temp = u_interval_order;
                        //idx = projectPoint.NbPoints();
                        flag = 1;
                    }



                }
                if (flag != 0)
                {
                    //cout << "最终挑选了： " << idx_max << endl<<endl<<endl;
                    mesh2surface_temp[k].energy = e_temp;
                    m2s_dist.push_back(e_temp);
                    m2s_grad_temp.push_back(grad_temp);
                    m2s_hessian_temp.push_back(B_temp);
                    mesh2surface_temp[k].u_basis_temp = basis_temp;
                    mesh2surface_temp[k].u_interval_order = u_interval_order_temp;
                    mesh2surface_temp[k].u = u_temp;
                    mesh2surface_temp[k].v = v_temp;
                    ioszz1 << mesh2surface_temp[k].point_on_mesh[0] << " " << mesh2surface_temp[k].point_on_mesh[1] << " " << mesh2surface_temp[k].point_on_mesh[2] << " " << bspline_surface_value1[0].value().value() << " " << bspline_surface_value1[1].value().value() << " " << bspline_surface_value1[2].value().value() << endl;
                    //cout << u_temp << endl;

                }
                else
                {
                    false_num++;
                    mesh2surface_temp[k].is_found = false;
                    cout << k << " " << mesh2surface_temp[k].u << endl;
                    continue;


                }
            }
            else
            {
                mesh2surface_temp[k].is_found = false;
                continue;
            }
        }
        else
        {
            double dot_temp = 0.0;
            int idx_max = 0;


            int flag = 0;
            double e_temp = 0;
            VectorXd grad_temp;

            Eigen::MatrixXd B_temp;
            B_temp.resize(ctr_dim, ctr_dim);
            vector<double> basis_temp;
            int u_interval_order_temp;
            Standard_Real dist_min = 1e6;


            int u_interval_order = u_k_position(mesh2surface_temp[k].u, u_knot_temp);


            vector<double> u_basis_temp;
            if (u_interval_order == 0)
            {
                vector<double> basis_temp(poly_degree + 1, 0);
                basis_temp[0] = 1;
                u_basis_temp = basis_temp;
            }
            else if (u_interval_order == u_knot_temp.size() - 1)
            {
                vector<double> basis_temp(poly_degree + 1, 0);
                basis_temp[poly_degree] = 1;
                u_basis_temp = basis_temp;
            }
            else
            {
                calc_basis_fuction(mesh2surface_temp[k].u, u_interval_order, poly_degree, u_basis_temp);
            }

            // cout << mesh2surface_temp[k].u << "  " << u_interval_order << endl;
            mesh2surface_temp[k].u_basis_temp = u_basis_temp;
            mesh2surface_temp[k].u_interval_order = u_interval_order;

            Vec3 bspline_surface_value;
            bspline_surface_value.setZero();

            if (u_interval_order == 0)
            {
                Vec3 pi0, pi1;
                pi0 = P[0][0];
                pi1 = P[0][1];
                bspline_surface_value[0] = mesh2surface_temp[k].v * pi1[0] + (1 - mesh2surface_temp[k].v) * pi0[0];
                bspline_surface_value[1] = mesh2surface_temp[k].v * pi1[1] + (1 - mesh2surface_temp[k].v) * pi0[1];
                bspline_surface_value[2] = mesh2surface_temp[k].v * pi1[2] + (1 - mesh2surface_temp[k].v) * pi0[2];
            }
            else if (u_interval_order == u_knot_temp.size() - 1)
            {
                Vec3 pi0, pi1;
                pi0 = P[curve_ctr_num - 1][0];
                pi1 = P[curve_ctr_num - 1][1];
                bspline_surface_value[0] = mesh2surface_temp[k].v * pi1[0] + (1 - mesh2surface_temp[k].v) * pi0[0];
                bspline_surface_value[1] = mesh2surface_temp[k].v * pi1[1] + (1 - mesh2surface_temp[k].v) * pi0[1];
                bspline_surface_value[2] = mesh2surface_temp[k].v * pi1[2] + (1 - mesh2surface_temp[k].v) * pi0[2];
            }
            else
            {
                for (int i = 0; i < u_basis_temp.size(); i++)
                {
                    Vec3 pi0, pi1;
                    pi0 = P[u_interval_order - poly_degree + i][0];
                    pi1 = P[u_interval_order - poly_degree + i][1];
                    bspline_surface_value[0] += u_basis_temp[i] * mesh2surface_temp[k].v * pi1[0] + u_basis_temp[i] * (1 - mesh2surface_temp[k].v) * pi0[0];
                    bspline_surface_value[1] += u_basis_temp[i] * mesh2surface_temp[k].v * pi1[1] + u_basis_temp[i] * (1 - mesh2surface_temp[k].v) * pi0[1];
                    bspline_surface_value[2] += u_basis_temp[i] * mesh2surface_temp[k].v * pi1[2] + u_basis_temp[i] * (1 - mesh2surface_temp[k].v) * pi0[2];
                }
            }

            scalar_t e = (-bspline_surface_value + mesh2surface_temp[k].point_on_mesh) * (-bspline_surface_value + mesh2surface_temp[k].point_on_mesh).transpose();

            //data_mesh2surface += e;



            Eigen::MatrixXd B;
            B.resize(ctr_dim, ctr_dim);

            for (int r = 0; r < ctr_dim; r++)
            {
                B.row(r) = e.derivatives()(r).derivatives().transpose();
            }
            B_temp = B;
            e_temp = e.value().value();
            grad_temp = e.value().derivatives();
            basis_temp = u_basis_temp;
            u_interval_order_temp = u_interval_order;
            flag = 1;





            m2s_dist.push_back(e_temp);
            m2s_grad_temp.push_back(grad_temp);
            m2s_hessian_temp.push_back(B_temp);
            mesh2surface_temp[k].u_basis_temp = basis_temp;
            mesh2surface_temp[k].u_interval_order = u_interval_order_temp;

        }
    }
    if (iter_num == 1)
    {
        double per = (double)false_num / (double)mesh2surface_temp.size();
        if (per > 0.15)
        {
            is_good = 0;
        }
    }
    cout <<"未找到投影点的采样点个数" << false_num << " " << mesh2surface_temp.size() << endl;
    if (iter_num == 0)
    {
        for (int k = 0; k < mesh2surface_temp.size(); k++)
        {


            if (nor_sum < 0)
            {
                mesh2surface_temp[k].normal = -mesh2surface_temp[k].normal;
                is_neg = -1;
            }

        }
    }
    double m2s_sum_dist = std::accumulate(m2s_dist.begin(), m2s_dist.end(), 0.0);
    double m2s_max = (*max_element(m2s_dist.begin(), m2s_dist.end()));


    double m2s_mean_dist = double(m2s_sum_dist) / double(m2s_dist.size());
    for (int d0 = 0; d0 < m2s_dist.size(); d0++)
    {
        double wt = std::pow(double(m2s_dist[d0]) / (double(m2s_mean_dist)), weight_order);
        if (iter_num >= 0)
        {
            wt = 1;
        }
        data_mesh2surface += wt * m2s_dist[d0];
        data_mesh2surface_grad += wt * m2s_grad_temp[d0];
        data_mesh2surface_hessian += wt * m2s_hessian_temp[d0];
        m2s_weight.push_back(wt);
    }
    /*
    scalar_t data_surface2mesh = 0;
    MatrixXd data_surface2mesh_hessian;
    VectorXd data_surface2mesh_grad;
    data_surface2mesh_grad.resize(ctr_dim);
    data_surface2mesh_grad.setZero();
    data_surface2mesh_hessian.resize(ctr_dim, ctr_dim);
    data_surface2mesh_hessian.setZero();

    for (int k = 0; k < surface2mesh_temp.size(); k++)
    {

        //auto it_u = std::upper_bound(u_knot_temp.begin(), u_knot_temp.end(), surface2mesh_temp[k].u);
        int u_interval_order = u_k_position(surface2mesh_temp[k].u,u_knot_temp);//std::distance(u_knot_temp.begin(), it_u) - 1;
        vector<double> u_basis_temp;

        if (u_interval_order == 0)
        {
            vector<double> basis_temp(poly_degree + 1, 0);
            basis_temp[0] = 1;
            u_basis_temp = basis_temp;
        }
        else if (u_interval_order == u_knot_temp.size() - 1)
        {
            vector<double> basis_temp(poly_degree + 1, 0);
            basis_temp[poly_degree] = 1;
            u_basis_temp = basis_temp;
        }
        else
        {
            calc_basis_fuction(mesh2surface_temp[k].u, u_interval_order, poly_degree, u_basis_temp);
        }

        Vec3 bspline_surface_value;
        bspline_surface_value.setZero();
        if (u_interval_order == 0)
        {
            Vec3 pi0, pi1;
            pi0 = P[0][0];
            pi1 = P[0][1];
            bspline_surface_value[0] = surface2mesh_temp[k].v * pi1[0] +  (1 - surface2mesh_temp[k].v) * pi0[0];
            bspline_surface_value[1] = surface2mesh_temp[k].v * pi1[1] +  (1 - surface2mesh_temp[k].v) * pi0[1];
            bspline_surface_value[2] = surface2mesh_temp[k].v * pi1[2] +  (1 - surface2mesh_temp[k].v) * pi0[2];
        }
        else if (u_interval_order == u_knot_temp.size()-1)
        {
            Vec3 pi0, pi1;
            pi0 = P[curve_ctr_num-1][0];
            pi1 = P[curve_ctr_num-1][1];
            bspline_surface_value[0] = surface2mesh_temp[k].v * pi1[0] + (1 - surface2mesh_temp[k].v) * pi0[0];
            bspline_surface_value[1] = surface2mesh_temp[k].v * pi1[1] + (1 - surface2mesh_temp[k].v) * pi0[1];
            bspline_surface_value[2] = surface2mesh_temp[k].v * pi1[2] + (1 - surface2mesh_temp[k].v) * pi0[2];
        }
        else
        {
            for (int i = 0; i < u_basis_temp.size(); i++)
            {
                Vec3 pi0, pi1;
                pi0 = P[u_interval_order - poly_degree + i][0];
                pi1 = P[u_interval_order - poly_degree + i][1];
                bspline_surface_value[0] += u_basis_temp[i] * surface2mesh_temp[k].v * pi1[0] + u_basis_temp[i] * (1 - surface2mesh_temp[k].v) * pi0[0];
                bspline_surface_value[1] += u_basis_temp[i] * surface2mesh_temp[k].v * pi1[1] + u_basis_temp[i] * (1 - surface2mesh_temp[k].v) * pi0[1];
                bspline_surface_value[2] += u_basis_temp[i] * surface2mesh_temp[k].v * pi1[2] + u_basis_temp[i] * (1 - surface2mesh_temp[k].v) * pi0[2];
            }
        }
        surface2mesh_temp[k].point_on_surface[0] = bspline_surface_value.x().value().value();
        surface2mesh_temp[k].point_on_surface[1] = bspline_surface_value.y().value().value();
        surface2mesh_temp[k].point_on_surface[2] = bspline_surface_value.z().value().value();

        auto result_closest = tree.closest_point_and_primitive(Point(surface2mesh_temp[k].point_on_surface[0].value().value(), surface2mesh_temp[k].point_on_surface[1].value().value(), surface2mesh_temp[k].point_on_surface[2].value().value()));
        surface2mesh_temp[k].point_on_mesh[0] = result_closest.first.x();
        surface2mesh_temp[k].point_on_mesh[1] = result_closest.first.y();
        surface2mesh_temp[k].point_on_mesh[2] = result_closest.first.z();

        scalar_t e = (bspline_surface_value - surface2mesh_temp[k].point_on_mesh) * (bspline_surface_value - surface2mesh_temp[k].point_on_mesh).transpose();
        //data_surface2mesh += e;

        s2m_dist.push_back(e.value().value());
        s2m_grad_temp.push_back(e.value().derivatives());

        //data_surface2mesh_grad += e.value().derivatives();


        Eigen::MatrixXd B;
        B.resize(ctr_dim, ctr_dim);
        for (int r = 0; r < ctr_dim; r++)
        {
            B.row(r) = e.derivatives()(r).derivatives().transpose();
        }

        s2m_hessian_temp.push_back(B);
        //data_surface2mesh_hessian += B;
    }


    double s2m_sum_dist = std::accumulate(s2m_dist.begin(), s2m_dist.end(), 0.0);

    double s2m_mean_dist = double(s2m_sum_dist) / double(s2m_dist.size());
    for (int d0 = 0; d0 < s2m_dist.size(); d0++)
    {
        double wt = std::pow(double(s2m_dist[d0]) / double(s2m_mean_dist), weight_order);
        data_surface2mesh += wt * s2m_dist[d0];
        data_surface2mesh_grad += wt * s2m_grad_temp[d0];
        data_surface2mesh_hessian += wt * s2m_hessian_temp[d0];
        s2m_weight.push_back(wt);
    }
    */
    data_value = data_mesh2surface_weight * data_mesh2surface;
    grad = data_mesh2surface_weight * data_mesh2surface_grad;
    hessian = data_mesh2surface_weight * data_mesh2surface_hessian;
    //cout << "梯度： " << data_surface2mesh_grad << endl;
    //cout << "海塞： " << hessian << endl;

}

void BSplineFitting::calc_smooth_term_grad_hessien(VectorXd& BSpline, scalar_t& smooth_value, Eigen::VectorXd& smooth_grad, Eigen::MatrixXd& smooth_hessian)
{
    smooth_value = 0;
    int num = BSpline.size();//X维数
    smooth_grad.resize(num);
    smooth_grad.setZero();
    smooth_hessian.resize(num, num);
    smooth_hessian.setZero();
    Vec12 X;
    X.resize(num);
    for (int r = 0; r < num / 3; r++)
    {
        X(3 * r).value() = BSpline(3 * r);
        X(3 * r + 1).value() = BSpline(3 * r + 1);
        X(3 * r + 2).value() = BSpline(3 * r + 2);
    }

    // 2.初始化一阶导数; repeat partial derivatives for the inner AutoDiffScalar
    for (int id = 0; id < num; id++)
    {
        X(id).derivatives().resize(num);
        X(id).derivatives().setZero();
        X(id).derivatives()(id) = 1;
        X(id).value().derivatives() = inner_derivative_t::Unit(num, id);
    }

    //  3.初始化海塞矩阵;set the hessian matrix to zero
    for (int idx = 0; idx < num; idx++)
    {
        for (int id = 0; id < num; id++)
        {
            X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(num);
        }
    }

    std::vector<vector<Vec3>> P;
    P.resize(num / 6);
    for (int i = 0; i < P.size(); i++)
    {
        P[i].resize(2);
    }

    for (int i = 0; i < P.size(); i++)
    {
        P[i][0][0] = X[3 * i];
        P[i][0][1] = X[3 * i + 1];
        P[i][0][2] = X[3 * i + 2];
        P[i][1][0] = X[3 * (num / 6) + 3 * i];
        P[i][1][1] = X[3 * (num / 6) + 3 * i + 1];
        P[i][1][2] = X[3 * (num / 6) + 3 * i + 2];
    }


    scalar_t s = 0;
    if (iter_num < 18)
    {
        for (int dim_col = 0; dim_col < 2; dim_col++)
        {
            for (int dim_row = 0; dim_row < num / 6 - 2; dim_row++)
            {
                Vec3 pij; pij = P[dim_row][dim_col];
                Vec3 pi1j; pi1j = P[dim_row + 1][dim_col];
                Vec3 pi2j; pi2j = P[dim_row + 2][dim_col];
                Vec3 ruls = (pi2j - pi1j) - (pi1j - pij);
                s += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
            }
        }

    }
    else
    {
        for (int dim = 0; dim < 2; dim++)
        {
            for (int k = 3; k < u_knot_temp.size() - 4; k++)
            {
                auto u = u_knot_temp;
                vector<vector<double>> B_const;
                B_const.resize(4);
                B_const[0].resize(1);
                B_const[1].resize(3);
                B_const[2].resize(3);
                B_const[3].resize(1);
                B_const[0][0] = 1.0 / ((u[k + 1] - u[k]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k - 2]));

                B_const[1][0] = 1.0 / ((u[k + 1] - u[k - 2]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[1][1] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[1][2] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[2][0] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[2][1] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[2][2] = 1.0 / ((u[k + 3] - u[k]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[3][0] = 1.0 / ((u[k + 3] - u[k]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                vector<vector<double>> sub_const;
                sub_const.resize(4);
                sub_const[0].resize(1);
                sub_const[1].resize(3);
                sub_const[2].resize(3);
                sub_const[3].resize(1);
                sub_const[0][0] = 6 * u[k + 1];

                sub_const[1][0] = -2 * u[k - 2] - 4 * u[k + 1];

                sub_const[1][1] = -2 * u[k + 2] - 2 * u[k + 1] - 2 * u[k - 1];

                sub_const[1][2] = -2 * u[k] - 4 * u[k + 2];

                sub_const[2][0] = 2 * u[k + 1] + 4 * u[k - 1];

                sub_const[2][1] = 2 * u[k - 1] + 2 * u[k] + 2 * u[k + 2];

                sub_const[2][2] = 4 * u[k] + 2 * u[k + 3];

                sub_const[3][0] = -6 * u[k];

                double delta = u[k + 1] - u[k];
                Vec3 p0 = P[k - 3][dim];
                Vec3 p1 = P[k - 2][dim];
                Vec3 p2 = P[k - 1][dim];
                Vec3 p3 = P[k][dim];
                s += (p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2]) * B_const[0][0] * B_const[0][0] * (12 * pow(delta, 3) - 3 * (sub_const[0][0] + sub_const[0][0]) * pow(delta, 2) + sub_const[0][0] * sub_const[0][0] * delta);

                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][0] * (12 * pow(delta, 3) + 3 * (sub_const[1][0] + sub_const[1][0]) * pow(delta, 2) + sub_const[1][0] * sub_const[1][0] * delta);
                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][1] * B_const[1][1] * (12 * pow(delta, 3) + 3 * (sub_const[1][1] + sub_const[1][1]) * pow(delta, 2) + sub_const[1][1] * sub_const[1][1] * delta);
                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][2] * B_const[1][2] * (12 * pow(delta, 3) + 3 * (sub_const[1][2] + sub_const[1][2]) * pow(delta, 2) + sub_const[1][2] * sub_const[1][2] * delta);
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][1] * (12 * pow(delta, 3) + 3 * (sub_const[1][0] + sub_const[1][1]) * pow(delta, 2) + (sub_const[1][0] * sub_const[1][1]) * delta);
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][2] * (12 * pow(delta, 3) + 3 * (sub_const[1][0] + sub_const[1][2]) * pow(delta, 2) + (sub_const[1][0] * sub_const[1][2]) * delta);
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][1] * B_const[1][2] * (12 * pow(delta, 3) + 3 * (sub_const[1][1] + sub_const[1][2]) * pow(delta, 2) + (sub_const[1][1] * sub_const[1][2]) * delta);

                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][0] * (12 * pow(delta, 3) - 3 * (sub_const[2][0] + sub_const[2][0]) * pow(delta, 2) + sub_const[2][0] * sub_const[2][0] * delta);
                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][1] * B_const[2][1] * (12 * pow(delta, 3) - 3 * (sub_const[2][1] + sub_const[2][1]) * pow(delta, 2) + sub_const[2][1] * sub_const[2][1] * delta);
                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][2] * B_const[2][2] * (12 * pow(delta, 3) - 3 * (sub_const[2][2] + sub_const[2][2]) * pow(delta, 2) + sub_const[2][2] * sub_const[2][2] * delta);
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][1] * (12 * pow(delta, 3) - 3 * (sub_const[2][0] + sub_const[2][1]) * pow(delta, 2) + (sub_const[2][0] * sub_const[2][1]) * delta);
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][2] * (12 * pow(delta, 3) - 3 * (sub_const[2][0] + sub_const[2][2]) * pow(delta, 2) + (sub_const[2][0] * sub_const[2][2]) * delta);
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][1] * B_const[2][2] * (12 * pow(delta, 3) - 3 * (sub_const[2][1] + sub_const[2][2]) * pow(delta, 2) + (sub_const[2][1] * sub_const[2][2]) * delta);

                s += (p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2]) * B_const[3][0] * B_const[3][0] * (12 * pow(delta, 3) + 3 * (sub_const[3][0] + sub_const[3][0]) * pow(delta, 2) + sub_const[3][0] * sub_const[3][0] * delta);

                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][0] * (-12 * pow(delta, 3) + 3 * (sub_const[0][0] - sub_const[1][0]) * pow(delta, 2) + (sub_const[0][0] * sub_const[1][0]) * delta);
                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][1] * (-12 * pow(delta, 3) + 3 * (sub_const[0][0] - sub_const[1][1]) * pow(delta, 2) + (sub_const[0][0] * sub_const[1][1]) * delta);
                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][2] * (-12 * pow(delta, 3) + 3 * (sub_const[0][0] - sub_const[1][2]) * pow(delta, 2) + (sub_const[0][0] * sub_const[1][2]) * delta);

                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][0] * (12 * pow(delta, 3) - 3 * (sub_const[0][0] + sub_const[2][0]) * pow(delta, 2) + (sub_const[0][0] * sub_const[2][0]) * delta);
                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][1] * (12 * pow(delta, 3) - 3 * (sub_const[0][0] + sub_const[2][1]) * pow(delta, 2) + (sub_const[0][0] * sub_const[2][1]) * delta);
                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][2] * (12 * pow(delta, 3) - 3 * (sub_const[0][0] + sub_const[2][2]) * pow(delta, 2) + (sub_const[0][0] * sub_const[2][2]) * delta);

                s += 2 * (p0[0] * p3[0] + p0[1] * p3[1] + p0[2] * p3[2]) * B_const[0][0] * B_const[3][0] * (-12 * pow(delta, 3) + 3 * (sub_const[0][0] - sub_const[3][0]) * pow(delta, 2) + (sub_const[0][0] * sub_const[3][0]) * delta);

                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][0] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][0] + sub_const[2][0]) * pow(delta, 2) + (sub_const[1][0] * sub_const[2][0]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][1] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][0] + sub_const[2][1]) * pow(delta, 2) + (sub_const[1][0] * sub_const[2][1]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][2] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][0] + sub_const[2][2]) * pow(delta, 2) + (sub_const[1][0] * sub_const[2][2]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][0] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][1] + sub_const[2][0]) * pow(delta, 2) + (sub_const[1][1] * sub_const[2][0]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][1] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][1] + sub_const[2][1]) * pow(delta, 2) + (sub_const[1][1] * sub_const[2][1]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][2] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][1] + sub_const[2][2]) * pow(delta, 2) + (sub_const[1][1] * sub_const[2][2]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][0] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][2] + sub_const[2][0]) * pow(delta, 2) + (sub_const[1][2] * sub_const[2][0]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][1] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][2] + sub_const[2][1]) * pow(delta, 2) + (sub_const[1][2] * sub_const[2][1]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][2] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][2] + sub_const[2][2]) * pow(delta, 2) + (sub_const[1][2] * sub_const[2][2]) * delta);

                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][0] * B_const[3][0] * (12 * pow(delta, 3) + 3 * (sub_const[1][0] + sub_const[3][0]) * pow(delta, 2) + (sub_const[1][0] * sub_const[3][0]) * delta);
                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][1] * B_const[3][0] * (12 * pow(delta, 3) + 3 * (sub_const[1][1] + sub_const[3][0]) * pow(delta, 2) + (sub_const[1][1] * sub_const[3][0]) * delta);
                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][2] * B_const[3][0] * (12 * pow(delta, 3) + 3 * (sub_const[1][2] + sub_const[3][0]) * pow(delta, 2) + (sub_const[1][2] * sub_const[3][0]) * delta);

                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][0] * B_const[3][0] * (-12 * pow(delta, 3) + 3 * (sub_const[2][0] - sub_const[3][0]) * pow(delta, 2) + (sub_const[2][0] * sub_const[3][0]) * delta);
                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][1] * B_const[3][0] * (-12 * pow(delta, 3) + 3 * (sub_const[2][1] - sub_const[3][0]) * pow(delta, 2) + (sub_const[2][1] * sub_const[3][0]) * delta);
                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][2] * B_const[3][0] * (-12 * pow(delta, 3) + 3 * (sub_const[2][2] - sub_const[3][0]) * pow(delta, 2) + (sub_const[2][2] * sub_const[3][0]) * delta);
            }
        }
        for (int dim = 0; dim < 0; dim++)
        {
            for (int k = 3; k < u_knot_temp.size() - 4; k++)
            {
                auto u = u_knot_temp;
                vector<vector<double>> B_const;
                B_const.resize(4);
                B_const[0].resize(1);
                B_const[1].resize(3);
                B_const[2].resize(3);
                B_const[3].resize(1);
                B_const[0][0] = 1.0 / ((u[k + 1] - u[k]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k - 2]));

                B_const[1][0] = 1.0 / ((u[k + 1] - u[k - 2]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[1][1] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[1][2] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[2][0] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[2][1] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[2][2] = 1.0 / ((u[k + 3] - u[k]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[3][0] = 1.0 / ((u[k + 3] - u[k]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));



                double delta = u[k + 1] - u[k];
                Vec3 p0 = P[k - 3][dim];
                Vec3 p1 = P[k - 2][dim];
                Vec3 p2 = P[k - 1][dim];
                Vec3 p3 = P[k][dim];
                s += (p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2]) * B_const[0][0] * B_const[0][0] * 36 * delta;

                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][0] * 36 * delta;
                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][1] * B_const[1][1] * 36 * delta;
                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][2] * B_const[1][2] * 36 * delta;
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][1] * 36 * delta;
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][2] * 36 * delta;
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][1] * B_const[1][2] * 36 * delta;

                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][0] * 36 * delta;
                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][1] * B_const[2][1] * 36 * delta;
                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][2] * B_const[2][2] * 36 * delta;
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][1] * 36 * delta;
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][2] * 36 * delta;
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][1] * B_const[2][2] * 36 * delta;

                s += (p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2]) * B_const[3][0] * B_const[3][0] * 36 * delta;

                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][0] * (-36) * delta;
                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][1] * (-36) * delta;
                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][2] * (-36) * delta;

                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][0] * 36 * delta;
                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][1] * 36 * delta;
                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][2] * 36 * delta;

                s += 2 * (p0[0] * p3[0] + p0[1] * p3[1] + p0[2] * p3[2]) * B_const[0][0] * B_const[3][0] * (-36) * delta;

                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][0] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][1] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][2] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][0] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][1] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][2] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][0] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][1] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][2] * (-36) * delta;

                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][0] * B_const[3][0] * 36 * delta;
                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][1] * B_const[3][0] * 36 * delta;
                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][2] * B_const[3][0] * 36 * delta;

                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][0] * B_const[3][0] * (-36) * delta;
                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][1] * B_const[3][0] * (-36) * delta;
                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][2] * B_const[3][0] * (-36) * delta;
            }
        }
        //s *= 1e-4;
    }

    scalar_t s1 = 0;

    for (int dim_row = 0; dim_row < num / 6 - 1; dim_row++)
    {
        Vec3 pi0; pi0 = P[dim_row][0];
        Vec3 pi1; pi1 = P[dim_row][1];

        Vec3 pii0; pii0 = P[dim_row + 1][0];
        Vec3 pii1; pii1 = P[dim_row + 1][1];

        Vec3 ruls = (pii1 - pii0) - (pi1 - pi0);
        s1 += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
    }

    vector<double> res = uni_para(5);
    vector<vector<Vec3>> curve_res;
    curve_res.resize(2);

    for (int i = 0; i < res.size(); i++)
    {
        int k = u_k_position(res[i], u_knot_temp);
        vector<double> u_basis_temp;
        calc_basis_fuction(res[i], k, poly_degree, u_basis_temp);
        Vec3 bspline_surface_value1;
        bspline_surface_value1.setZero();
        Vec3 bspline_surface_value2;
        bspline_surface_value2.setZero();
        if (k == 0)
        {

            bspline_surface_value1 = P[0][0];

            bspline_surface_value2 = P[0][1];
        }
        else if (k == u_knot_temp.size() - 1)
        {
            bspline_surface_value1 = P[BSpline.size() / 6 - 1][0];

            bspline_surface_value2 = P[BSpline.size() / 6 - 1][1];
        }
        else
        {
            for (int j = 0; j < u_basis_temp.size(); j++)
            {
                bspline_surface_value1[0] += u_basis_temp[j] * P[k - 3 + j][0][0];
                bspline_surface_value1[1] += u_basis_temp[j] * P[k - 3 + j][0][1];
                bspline_surface_value1[2] += u_basis_temp[j] * P[k - 3 + j][0][2];
                bspline_surface_value2[0] += u_basis_temp[j] * P[k - 3 + j][1][0];
                bspline_surface_value2[1] += u_basis_temp[j] * P[k - 3 + j][1][1];
                bspline_surface_value2[2] += u_basis_temp[j] * P[k - 3 + j][1][2];
            }
        }
        curve_res[0].push_back(bspline_surface_value1);
        curve_res[1].push_back(bspline_surface_value2);

    }

    for (int dim_row = 0; dim_row < curve_res[0].size() - 3; dim_row++)
    {
        Vec3 pi0 = curve_res[0][dim_row];
        Vec3 pi1 = curve_res[1][dim_row];

        Vec3 pii0; pii0 = curve_res[0][dim_row + 3];
        Vec3 pii1; pii1 = curve_res[1][dim_row + 3];

        Vec3 ruls = (pii1 - pii0) - (pi1 - pi0);
        s1 += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
    }
    for (int dim_col = 0; dim_col < 2; dim_col++)
    {
        for (int dim_row = 0; dim_row < curve_res[0].size() - 4; dim_row++)
        {
            Vec3 pi0 = curve_res[dim_col][dim_row];
            Vec3 pi1 = curve_res[dim_col][dim_row + 2];

            Vec3 pii0; pii0 = curve_res[dim_col][dim_row + 4];
            Vec3 ruls = pi0 - pi1 + pii0 - pi1;
            s += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
        }
    }

    scalar_t s_total = u_smooth_weight * s + v_smooth_weight * s1;
    smooth_value = s_total;
    smooth_grad = s_total.value().derivatives();
    Eigen::MatrixXd B;
    B.resize(num, num);
    for (int r = 0; r < num; r++)
    {
        B.row(r) = s_total.derivatives()(r).derivatives().transpose();
    }
    smooth_hessian += B;
}

void BSplineFitting::calc_objective_grad_hessien(const Data& spline, scalar_t& data_value, Eigen::VectorXd& data_grad, Eigen::MatrixXd& data_hessian, Eigen::VectorXd& decrease_direction, scalar_t& smooth_value, Eigen::VectorXd& smooth_grad, Eigen::MatrixXd& smooth_hessian, Eigen::VectorXd& total_grad, Eigen::MatrixXd& total_hessian)
{

    VectorXd bspline;
    bspline.resize(spline.rows() * 3);
    for (int i = 0; i < spline.rows(); i++)
    {

        bspline[3 * i] = spline.coeff(i, 0);
        bspline[3 * i + 1] = spline.coeff(i, 1);
        bspline[3 * i + 2] = spline.coeff(i, 2);
    }

    smooth_value = 0;
    int ctr_dim = bspline.size();//X维数
    smooth_grad.resize(ctr_dim);
    smooth_grad.setZero();
    smooth_hessian.resize(ctr_dim, ctr_dim);
    smooth_hessian.setZero();

    decrease_direction.resize(ctr_dim);
    decrease_direction.setZero();
    total_grad.resize(ctr_dim);
    total_grad.setZero();
    total_hessian.resize(ctr_dim, ctr_dim);
    total_hessian.setZero();

    calc_data_term_grad_hessien(bspline, data_value, data_grad, data_hessian);

    calc_smooth_term_grad_hessien(bspline, smooth_value, smooth_grad, smooth_hessian);
    cout << "data" << data_grad.norm() << endl;
    cout << "smooth" << smooth_grad.norm() << endl;



    //total_grad = data_grad;
    //total_hessian = data_hessian;

    // cout << data_grad<<endl;
    //calc_smooth_term_grad_hessien(bspline, smooth_value, smooth_grad, smooth_hessian);
    //total_grad = data_grad + smooth_grad;
    //total_hessian = data_hessian + smooth_hessian;

    /*
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(total_hessian, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::VectorXd singularValues = svd.singularValues();
    Eigen::VectorXd inverse_singularValues;
    inverse_singularValues.resize(singularValues.size());
    inverse_singularValues.setZero();
    for (int i = 0; i < singularValues.size(); i++)//4.62223e-33数值不稳定
    {
        if (singularValues[i] < 1e-6)
        {
            singularValues[i] = 1e-6;
            inverse_singularValues[i] = 1.0 / singularValues[i];
        }
        else
        {
            inverse_singularValues[i] = 1.0 / singularValues[i];
        }
    }

    auto pseudo_inverse_matrix = V * inverse_singularValues.asDiagonal() * U.transpose();
    decrease_direction = -pseudo_inverse_matrix * total_grad;
    decrease_direction.normalize();
    */
}

double BSplineFitting::calc_data_term_energy(const Data& spline)
{
    VectorXd bspline;
    bspline.resize(spline.rows() * 3);
    for (int i = 0; i < spline.rows(); i++)
    {

        bspline[3 * i] = spline.coeff(i, 0);
        bspline[3 * i + 1] = spline.coeff(i, 1);
        bspline[3 * i + 2] = spline.coeff(i, 2);
    }

    double data_value = 0.0;
    int ctr_dim = bspline.size();
    int curve_ctr_num = ctr_dim / 6;
    std::vector<vector<Vec3>> P;
    P.resize(curve_ctr_num);
    for (int i = 0; i < P.size(); i++)
    {
        P[i].resize(2);
    }

    for (int i = 0; i < P.size(); i++)
    {
        P[i][0][0] = bspline[3 * i];
        P[i][0][1] = bspline[3 * i + 1];
        P[i][0][2] = bspline[3 * i + 2];
        P[i][1][0] = bspline[3 * curve_ctr_num + 3 * i];
        P[i][1][1] = bspline[3 * curve_ctr_num + 3 * i + 1];
        P[i][1][2] = bspline[3 * curve_ctr_num + 3 * i + 2];
    }

    double data_mesh2surface_value = 0.0;

    vector<CtrInfo> ctr_point1;
    for (int j = 0; j < curve_ctr_num; j++)
    {
        CtrInfo ci;
        Standard_Real px = bspline[3 * j];
        Standard_Real py = bspline[3 * j + 1];
        Standard_Real pz = bspline[3 * j + 2];
        ci.position = { px,py,pz };
        ci.weight = 1;
        ctr_point1.push_back(ci);
    }

    vector<CtrInfo> ctr_point2;
    for (int j = 0; j < curve_ctr_num; j++)
    {
        CtrInfo ci;
        Standard_Real px = bspline[3 * (j + curve_ctr_num)];
        Standard_Real py = bspline[3 * (j + curve_ctr_num) + 1];
        Standard_Real pz = bspline[3 * (j + curve_ctr_num) + 2];
        ci.position = { px,py,pz };
        ci.weight = 1;
        ctr_point2.push_back(ci);
    }

    Handle(Geom_BSplineSurface) BSpline_surface = nullptr;
    create_BSpline_surface(ctr_point1, ctr_point2, BSpline_surface);


    vector<double> m2s_energy; m2s_energy.clear();
    m2s_energy.shrink_to_fit();
    //iter_num++;
    //BSpline_surface_viewer_2(spline, 100, iter_num);
    //ofstream ioszz("mesh2surface" + to_string(iter_num) + ".txt");
    e_intervals.resize(u_knot_temp.size());
    avr_e_intervals.resize(u_knot_temp.size());
    e_len.resize(u_knot_temp.size());
    e_num.resize(u_knot_temp.size());
    for (int i = 0; i < e_intervals.size(); i++)
    {
        e_intervals[i] = 0;
        avr_e_intervals[i] = 0;
        e_num[i] = 0;
        e_len[i] = 0;
    }

    max_energy = 0.0;
    min_energy = 10.0;
    for (int k = 0; k < mesh2surface_temp.size(); k++)
    {
        if (mesh2surface_temp[k].is_found)
        {
            if (((change_idx - mesh2surface_temp[k].u_interval_order) <= 0) && ((change_idx - mesh2surface_temp[k].u_interval_order) > -4))
            {
                Vec3 bspline_surface_value;
                bspline_surface_value.setZero();
                vector<double> u_basis_temp = mesh2surface_temp[k].u_basis_temp;
                int u_interval_order = mesh2surface_temp[k].u_interval_order;
                if (u_interval_order == 0)
                {
                    Vec3 pi0, pi1;
                    pi0 = P[0][0];
                    pi1 = P[0][1];
                    bspline_surface_value[0] = mesh2surface_temp[k].v * pi1[0] + (1 - mesh2surface_temp[k].v) * pi0[0];
                    bspline_surface_value[1] = mesh2surface_temp[k].v * pi1[1] + (1 - mesh2surface_temp[k].v) * pi0[1];
                    bspline_surface_value[2] = mesh2surface_temp[k].v * pi1[2] + (1 - mesh2surface_temp[k].v) * pi0[2];
                }
                else if (u_interval_order == u_knot_temp.size() - 1)
                {
                    Vec3 pi0, pi1;
                    pi0 = P[curve_ctr_num - 1][0];
                    pi1 = P[curve_ctr_num - 1][1];
                    bspline_surface_value[0] = mesh2surface_temp[k].v * pi1[0] + (1 - mesh2surface_temp[k].v) * pi0[0];
                    bspline_surface_value[1] = mesh2surface_temp[k].v * pi1[1] + (1 - mesh2surface_temp[k].v) * pi0[1];
                    bspline_surface_value[2] = mesh2surface_temp[k].v * pi1[2] + (1 - mesh2surface_temp[k].v) * pi0[2];
                }
                else
                {
                    for (int i = 0; i < u_basis_temp.size(); i++)
                    {
                        Vec3 pi0, pi1;
                        pi0 = P[u_interval_order - poly_degree + i][0];
                        pi1 = P[u_interval_order - poly_degree + i][1];
                        bspline_surface_value[0] += u_basis_temp[i] * mesh2surface_temp[k].v * pi1[0] + u_basis_temp[i] * (1 - mesh2surface_temp[k].v) * pi0[0];
                        bspline_surface_value[1] += u_basis_temp[i] * mesh2surface_temp[k].v * pi1[1] + u_basis_temp[i] * (1 - mesh2surface_temp[k].v) * pi0[1];
                        bspline_surface_value[2] += u_basis_temp[i] * mesh2surface_temp[k].v * pi1[2] + u_basis_temp[i] * (1 - mesh2surface_temp[k].v) * pi0[2];
                    }
                }

                //ioszz << bspline_surface_value[0].value().value() << " " << bspline_surface_value[1].value().value() << " " << bspline_surface_value[2].value().value() << " " << mesh2surface_temp[k].point_on_mesh[0].value().value() << " " << mesh2surface_temp[k].point_on_mesh[1].value().value() << " " << mesh2surface_temp[k].point_on_mesh[2].value().value() << endl;
                //gp_Pnt tpp;
                //BSpline_surface->D0(mesh2surface_temp[k].u, mesh2surface_temp[k].v, tpp);

                //cout << tpp.X() << " " << tpp.Y() << " " << tpp.Z() << endl;
                //cout << bspline_surface_value.x().value().value() << " " << bspline_surface_value.y().value().value() << " " << bspline_surface_value.z().value().value() << endl;
                scalar_t e = (bspline_surface_value - mesh2surface_temp[k].point_on_mesh) * (bspline_surface_value - mesh2surface_temp[k].point_on_mesh).transpose();
                mesh2surface_temp[k].energy = e.value().value();
            }
            m2s_energy.push_back(mesh2surface_temp[k].energy);

            //cout << mesh2surface_temp[k].u << " " << mesh2surface_temp[k].u_interval_order << endl;
            e_intervals[mesh2surface_temp[k].u_interval_order] += mesh2surface_temp[k].energy;
            if (mesh2surface_temp[k].energy > max_energy)
            {
                max_energy = mesh2surface_temp[k].energy;
            }
            if (mesh2surface_temp[k].energy < min_energy)
            {
                min_energy = mesh2surface_temp[k].energy;
            }
            e_num[mesh2surface_temp[k].u_interval_order]++;
            //data_mesh2surface_value += e.value().value();
        }
    }
    for (int i = 0; i < e_intervals.size(); i++)
    {
        if (e_num[i] != 0)
        {
            avr_e_intervals[i] = e_intervals[i]/e_num[i];
        }
    }
    insert1_idx = -1;
    double max_temp_e = 0;
    for (int i = 0; i < e_intervals.size() - 1; i++)
    {
        double loc1 = u_knot_temp[i];
        double loc2 = u_knot_temp[i + 1];
        int u_interval_order = u_k_position(loc1, u_knot_temp);


        vector<double> u_basis_temp;
        if (u_interval_order == 0)
        {
            vector<double> basis_temp(poly_degree + 1, 0);
            basis_temp[0] = 1;
            u_basis_temp = basis_temp;
        }
        else if (u_interval_order == u_knot_temp.size() - 1)
        {
            vector<double> basis_temp(poly_degree + 1, 0);
            basis_temp[poly_degree] = 1;
            u_basis_temp = basis_temp;
        }
        else
        {
            calc_basis_fuction(loc1, u_interval_order, poly_degree, u_basis_temp);
        }

        Vec3 bspline_surface_value;
        bspline_surface_value.setZero();

        if (u_interval_order == 0)
        {
            Vec3 pi0, pi1;
            pi0 = P[0][0];
            pi1 = P[0][1];
            bspline_surface_value[0] = 0.5 * pi1[0] + 0.5 * pi0[0];
            bspline_surface_value[1] = 0.5 * pi1[1] + 0.5 * pi0[1];
            bspline_surface_value[2] = 0.5 * pi1[2] + 0.5 * pi0[2];
        }
        else if (u_interval_order == u_knot_temp.size() - 1)
        {
            Vec3 pi0, pi1;
            pi0 = P[curve_ctr_num - 1][0];
            pi1 = P[curve_ctr_num - 1][1];
            bspline_surface_value[0] = 0.5 * pi1[0] + 0.5 * pi0[0];
            bspline_surface_value[1] = 0.5 * pi1[1] + 0.5 * pi0[1];
            bspline_surface_value[2] = 0.5 * pi1[2] + 0.5 * pi0[2];
        }
        else
        {
            for (int i = 0; i < u_basis_temp.size(); i++)
            {
                Vec3 pi0, pi1;
                pi0 = P[u_interval_order - poly_degree + i][0];
                pi1 = P[u_interval_order - poly_degree + i][1];
                bspline_surface_value[0] += u_basis_temp[i] * 0.5 * pi1[0] + u_basis_temp[i] * 0.5 * pi0[0];
                bspline_surface_value[1] += u_basis_temp[i] * 0.5 * pi1[1] + u_basis_temp[i] * 0.5 * pi0[1];
                bspline_surface_value[2] += u_basis_temp[i] * 0.5 * pi1[2] + u_basis_temp[i] * 0.5 * pi0[2];
            }
        }
        int u_interval_order1 = u_k_position(loc2, u_knot_temp);
        vector<double> u_basis_temp1;
        if (u_interval_order1 == 0)
        {
            vector<double> basis_temp1(poly_degree + 1, 0);
            basis_temp1[0] = 1;
            u_basis_temp1 = basis_temp1;
        }
        else if (u_interval_order1 == u_knot_temp.size() - 1)
        {
            vector<double> basis_temp1(poly_degree + 1, 0);
            basis_temp1[poly_degree] = 1;
            u_basis_temp1 = basis_temp1;
        }
        else
        {
            calc_basis_fuction(loc2, u_interval_order1, poly_degree, u_basis_temp1);
        }

        Vec3 bspline_surface_value1;
        bspline_surface_value1.setZero();

        if (u_interval_order1 == 0)
        {
            Vec3 pi0, pi1;
            pi0 = P[0][0];
            pi1 = P[0][1];
            bspline_surface_value1[0] = 0.5 * pi1[0] + 0.5 * pi0[0];
            bspline_surface_value1[1] = 0.5 * pi1[1] + 0.5 * pi0[1];
            bspline_surface_value1[2] = 0.5 * pi1[2] + 0.5 * pi0[2];
        }
        else if (u_interval_order1 == u_knot_temp.size() - 1)
        {
            Vec3 pi0, pi1;
            pi0 = P[curve_ctr_num - 1][0];
            pi1 = P[curve_ctr_num - 1][1];
            bspline_surface_value1[0] = 0.5 * pi1[0] + 0.5 * pi0[0];
            bspline_surface_value1[1] = 0.5 * pi1[1] + 0.5 * pi0[1];
            bspline_surface_value1[2] = 0.5 * pi1[2] + 0.5 * pi0[2];
        }
        else
        {
            for (int i = 0; i < u_basis_temp1.size(); i++)
            {
                Vec3 pi0, pi1;
                pi0 = P[u_interval_order1 - poly_degree + i][0];
                pi1 = P[u_interval_order1 - poly_degree + i][1];
                bspline_surface_value1[0] += u_basis_temp1[i] * 0.5 * pi1[0] + u_basis_temp1[i] * 0.5 * pi0[0];
                bspline_surface_value1[1] += u_basis_temp1[i] * 0.5 * pi1[1] + u_basis_temp1[i] * 0.5 * pi0[1];
                bspline_surface_value1[2] += u_basis_temp1[i] * 0.5 * pi1[2] + u_basis_temp1[i] * 0.5 * pi0[2];
            }
        }
        scalar_t int_leng = (bspline_surface_value - bspline_surface_value1) * (bspline_surface_value - bspline_surface_value1).transpose();
        double leng1 = int_leng.value().value();
        leng1 = sqrt(leng1);
        cout << "i= " << i << "; " << "[" << u_knot_temp[i] << "," << u_knot_temp[i + 1] << "]" << e_intervals[i]<<" "<<avr_e_intervals[i] << ";length="<<leng1  << endl;
        e_len[i] = leng1;
       
        int init_interval = floor(u_knot_temp[i]);
        if (max_temp_e < e_intervals[i] && count_divide[init_interval] < 4&&leng1>0.02 && (count_divide[init_interval] < 1 || e_intervals[i]>0.1))
        {

            max_temp_e = e_intervals[i];
            if (u_knot_temp[i + 1] == 0)
            {
                int i_ = i;
                while (u_knot_temp[i_ + 1] == 0)
                {
                    i_++;
                }
                insert1_idx = i_;
            }
            else if (u_knot_temp[i] == u_knot_temp.back())
            {
                int i_ = i;
                while (u_knot_temp[i_] == u_knot_temp.back())
                {
                    i_--;
                }
                insert1_idx = i_;
            }
            else
            {
                insert1_idx = i;
            }

        }
    }


    if (max_temp_e == 0)
    {
        for (int i = 0; i < e_intervals.size() - 1; i++)
        {
            int init_interval = floor(u_knot_temp[i]);




            if (max_temp_e < e_intervals[i]&&e_len[i]>0.02 && (count_divide[init_interval] < 1 || avr_e_intervals[i]>4e-4))
            {

                max_temp_e = e_intervals[i];
                if (u_knot_temp[i + 1] == 0)
                {
                    int i_ = i;
                    while (u_knot_temp[i_ + 1] == 0)
                    {
                        i_++;
                    }
                    insert1_idx = i_;
                }
                else if (u_knot_temp[i] == u_knot_temp.back())
                {
                    int i_ = i;
                    while (u_knot_temp[i_] == u_knot_temp.back())
                    {
                        i_--;
                    }
                    insert1_idx = i_;
                }
                else
                {
                    insert1_idx = i;
                }

            }
        }
    }

    if (insert1_idx != -1)
    {
        cout << "插入区间：" << insert1_idx << " " << "[" << u_knot_temp[insert1_idx] << "," << u_knot_temp[insert1_idx + 1] << "]" << endl;


        cout << "区间细分次数" << count_divide[(int)floor(u_knot_temp[insert1_idx])] << endl;
    }
    else
    {
        cout << "不插入新控制点" << endl;
    }
    //double m2s_sum_dist = std::accumulate(m2s_energy.begin(), m2s_energy.end(), 0.0);
   //double m2s_mean_dist = double(m2s_sum_dist) / double(m2s_energy.size());
    for (int d0 = 0; d0 < m2s_energy.size(); d0++)
    {
        //double wt = std::pow(double(m2s_energy[d0]) / double(m2s_mean_dist), weight_order);
        data_mesh2surface_value += m2s_energy[d0];

    }
    avr_energy = data_mesh2surface_value / m2s_energy.size();

    //ioszz.close();

    double data_surface2mesh = 0.0;


    //iter_num++;
    //BSpline_surface_viewer_2(spline, 100, iter_num);
    //ofstream ioszz("surface2mesh" + to_string(iter_num) + ".txt");
    /*
    for (int k = 0; k < surface2mesh_temp.size(); k++)
    {
        gp_Pnt P_on_surface;
        BSpline_surface->D0(surface2mesh_temp[k].u, surface2mesh_temp[k].v, P_on_surface);
        Vector3d surface_point = { P_on_surface.X(), P_on_surface.Y(), P_on_surface.Z() };
        auto result_closest = tree.closest_point_and_primitive(Point(surface_point[0], surface_point[1], surface_point[2]));
        Vector3d close_point;
        close_point[0] = result_closest.first.x();
        close_point[1] = result_closest.first.y();
        close_point[2] = result_closest.first.z();

        double e = (surface_point - close_point).transpose() * (surface_point - close_point);
        //data_surface2mesh += e;
        s2m_energy.push_back(e);
        //ioszz << surface_point[0] << " " << surface_point[1] << " " << surface_point[2] << " " << close_point[0] << " " << close_point[1] << " " << close_point[2] << endl;
    }
    //ioszz.close();

    //double s2m_sum_dist = std::accumulate(s2m_energy.begin(), s2m_energy.end(), 0.0);
    //double s2m_mean_dist = double(s2m_sum_dist) / double(s2m_energy.size());
    for (int d0 = 0; d0 < s2m_energy.size(); d0++)
    {
        //double wt = std::pow(double(s2m_energy[d0]) / double(s2m_mean_dist), weight_order);
        data_surface2mesh += s2m_weight[d0] *  s2m_energy[d0];
    }
    */
    data_value = data_mesh2surface_weight * data_mesh2surface_value + data_surface2mesh_weight * data_surface2mesh;
    std::cout << "data_mesh2surface_weight * data_mesh2surface_value = " << data_mesh2surface_weight * data_mesh2surface_value << endl;
    std::cout << "data_surface2mesh_weight * data_surface2mesh = " << data_surface2mesh_weight * data_surface2mesh << endl << endl;;
    return data_value;
}

double BSplineFitting::calc_smooth_term_energyn(const Data& spline)
{
    VectorXd spline_v;
    spline_v.resize(spline.rows() * 3);
    int curve_order_num = (spline.rows() / 2) - 1;
    for (int i = 0; i < spline.rows(); i++)
    {
        spline_v[3 * i] = spline.coeff(i, 0);
        spline_v[3 * i + 1] = spline.coeff(i, 1);
        spline_v[3 * i + 2] = spline.coeff(i, 2);
    }

    double data_value = 0.0;
    std::vector<vector<Vec3>> P;
    P.resize((curve_order_num + 1));
    for (int i = 0; i < P.size(); i++)
    {
        P[i].resize(2);
    }

    for (int i = 0; i < P.size(); i++)
    {
        P[i][0][0] = spline_v[3 * i];
        P[i][0][1] = spline_v[3 * i + 1];
        P[i][0][2] = spline_v[3 * i + 2];
        P[i][1][0] = spline_v[3 * (curve_order_num + 1) + 3 * i];
        P[i][1][1] = spline_v[3 * (curve_order_num + 1) + 3 * i + 1];
        P[i][1][2] = spline_v[3 * (curve_order_num + 1) + 3 * i + 2];
    }

    double smooth_value = 0.0;
    scalar_t s = 0;
    if (iter_num < 18 && flag123 == 0)
    {
        for (int dim_col = 0; dim_col < 2; dim_col++)
        {
            for (int dim_row = 0; dim_row < curve_order_num - 1; dim_row++)
            {
                Vec3 pij; pij = P[dim_row][dim_col];
                Vec3 pi1j; pi1j = P[dim_row + 1][dim_col];
                Vec3 pi2j; pi2j = P[dim_row + 2][dim_col];
                Vec3 ruls = (pi2j - pi1j) - (pi1j - pij);
                s += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
            }
        }
    }
    else
    {
        cout << "u_weight=" << u_smooth_weight << endl;
        for (int dim = 0; dim < 2; dim++)
        {
            for (int k = 3; k < u_knot_temp.size() - 4; k++)
            {
                auto u = u_knot_temp;
                vector<vector<double>> B_const;
                B_const.resize(4);
                B_const[0].resize(1);
                B_const[1].resize(3);
                B_const[2].resize(3);
                B_const[3].resize(1);
                B_const[0][0] = 1.0 / ((u[k + 1] - u[k]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k - 2]));

                B_const[1][0] = 1.0 / ((u[k + 1] - u[k - 2]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[1][1] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[1][2] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[2][0] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[2][1] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[2][2] = 1.0 / ((u[k + 3] - u[k]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[3][0] = 1.0 / ((u[k + 3] - u[k]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                vector<vector<double>> sub_const;
                sub_const.resize(4);
                sub_const[0].resize(1);
                sub_const[1].resize(3);
                sub_const[2].resize(3);
                sub_const[3].resize(1);
                sub_const[0][0] = 6 * u[k + 1];

                sub_const[1][0] = -2 * u[k - 2] - 4 * u[k + 1];

                sub_const[1][1] = -2 * u[k + 2] - 2 * u[k + 1] - 2 * u[k - 1];

                sub_const[1][2] = -2 * u[k] - 4 * u[k + 2];

                sub_const[2][0] = 2 * u[k + 1] + 4 * u[k - 1];

                sub_const[2][1] = 2 * u[k - 1] + 2 * u[k] + 2 * u[k + 2];

                sub_const[2][2] = 4 * u[k] + 2 * u[k + 3];

                sub_const[3][0] = -6 * u[k];

                double delta = u[k + 1] - u[k];
                Vec3 p0 = P[k - 3][dim];
                Vec3 p1 = P[k - 2][dim];
                Vec3 p2 = P[k - 1][dim];
                Vec3 p3 = P[k][dim];
                s += (p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2]) * B_const[0][0] * B_const[0][0] * (12 * pow(delta, 3) - 3 * (sub_const[0][0] + sub_const[0][0]) * pow(delta, 2) + sub_const[0][0] * sub_const[0][0] * delta);

                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][0] * (12 * pow(delta, 3) + 3 * (sub_const[1][0] + sub_const[1][0]) * pow(delta, 2) + sub_const[1][0] * sub_const[1][0] * delta);
                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][1] * B_const[1][1] * (12 * pow(delta, 3) + 3 * (sub_const[1][1] + sub_const[1][1]) * pow(delta, 2) + sub_const[1][1] * sub_const[1][1] * delta);
                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][2] * B_const[1][2] * (12 * pow(delta, 3) + 3 * (sub_const[1][2] + sub_const[1][2]) * pow(delta, 2) + sub_const[1][2] * sub_const[1][2] * delta);
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][1] * (12 * pow(delta, 3) + 3 * (sub_const[1][0] + sub_const[1][1]) * pow(delta, 2) + (sub_const[1][0] * sub_const[1][1]) * delta);
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][2] * (12 * pow(delta, 3) + 3 * (sub_const[1][0] + sub_const[1][2]) * pow(delta, 2) + (sub_const[1][0] * sub_const[1][2]) * delta);
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][1] * B_const[1][2] * (12 * pow(delta, 3) + 3 * (sub_const[1][1] + sub_const[1][2]) * pow(delta, 2) + (sub_const[1][1] * sub_const[1][2]) * delta);

                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][0] * (12 * pow(delta, 3) - 3 * (sub_const[2][0] + sub_const[2][0]) * pow(delta, 2) + sub_const[2][0] * sub_const[2][0] * delta);
                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][1] * B_const[2][1] * (12 * pow(delta, 3) - 3 * (sub_const[2][1] + sub_const[2][1]) * pow(delta, 2) + sub_const[2][1] * sub_const[2][1] * delta);
                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][2] * B_const[2][2] * (12 * pow(delta, 3) - 3 * (sub_const[2][2] + sub_const[2][2]) * pow(delta, 2) + sub_const[2][2] * sub_const[2][2] * delta);
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][1] * (12 * pow(delta, 3) - 3 * (sub_const[2][0] + sub_const[2][1]) * pow(delta, 2) + (sub_const[2][0] * sub_const[2][1]) * delta);
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][2] * (12 * pow(delta, 3) - 3 * (sub_const[2][0] + sub_const[2][2]) * pow(delta, 2) + (sub_const[2][0] * sub_const[2][2]) * delta);
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][1] * B_const[2][2] * (12 * pow(delta, 3) - 3 * (sub_const[2][1] + sub_const[2][2]) * pow(delta, 2) + (sub_const[2][1] * sub_const[2][2]) * delta);

                s += (p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2]) * B_const[3][0] * B_const[3][0] * (12 * pow(delta, 3) + 3 * (sub_const[3][0] + sub_const[3][0]) * pow(delta, 2) + sub_const[3][0] * sub_const[3][0] * delta);

                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][0] * (-12 * pow(delta, 3) + 3 * (sub_const[0][0] - sub_const[1][0]) * pow(delta, 2) + (sub_const[0][0] * sub_const[1][0]) * delta);
                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][1] * (-12 * pow(delta, 3) + 3 * (sub_const[0][0] - sub_const[1][1]) * pow(delta, 2) + (sub_const[0][0] * sub_const[1][1]) * delta);
                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][2] * (-12 * pow(delta, 3) + 3 * (sub_const[0][0] - sub_const[1][2]) * pow(delta, 2) + (sub_const[0][0] * sub_const[1][2]) * delta);

                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][0] * (12 * pow(delta, 3) - 3 * (sub_const[0][0] + sub_const[2][0]) * pow(delta, 2) + (sub_const[0][0] * sub_const[2][0]) * delta);
                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][1] * (12 * pow(delta, 3) - 3 * (sub_const[0][0] + sub_const[2][1]) * pow(delta, 2) + (sub_const[0][0] * sub_const[2][1]) * delta);
                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][2] * (12 * pow(delta, 3) - 3 * (sub_const[0][0] + sub_const[2][2]) * pow(delta, 2) + (sub_const[0][0] * sub_const[2][2]) * delta);

                s += 2 * (p0[0] * p3[0] + p0[1] * p3[1] + p0[2] * p3[2]) * B_const[0][0] * B_const[3][0] * (-12 * pow(delta, 3) + 3 * (sub_const[0][0] - sub_const[3][0]) * pow(delta, 2) + (sub_const[0][0] * sub_const[3][0]) * delta);

                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][0] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][0] + sub_const[2][0]) * pow(delta, 2) + (sub_const[1][0] * sub_const[2][0]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][1] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][0] + sub_const[2][1]) * pow(delta, 2) + (sub_const[1][0] * sub_const[2][1]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][2] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][0] + sub_const[2][2]) * pow(delta, 2) + (sub_const[1][0] * sub_const[2][2]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][0] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][1] + sub_const[2][0]) * pow(delta, 2) + (sub_const[1][1] * sub_const[2][0]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][1] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][1] + sub_const[2][1]) * pow(delta, 2) + (sub_const[1][1] * sub_const[2][1]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][2] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][1] + sub_const[2][2]) * pow(delta, 2) + (sub_const[1][1] * sub_const[2][2]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][0] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][2] + sub_const[2][0]) * pow(delta, 2) + (sub_const[1][2] * sub_const[2][0]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][1] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][2] + sub_const[2][1]) * pow(delta, 2) + (sub_const[1][2] * sub_const[2][1]) * delta);
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][2] * (-12 * pow(delta, 3) + 3 * (-sub_const[1][2] + sub_const[2][2]) * pow(delta, 2) + (sub_const[1][2] * sub_const[2][2]) * delta);

                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][0] * B_const[3][0] * (12 * pow(delta, 3) + 3 * (sub_const[1][0] + sub_const[3][0]) * pow(delta, 2) + (sub_const[1][0] * sub_const[3][0]) * delta);
                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][1] * B_const[3][0] * (12 * pow(delta, 3) + 3 * (sub_const[1][1] + sub_const[3][0]) * pow(delta, 2) + (sub_const[1][1] * sub_const[3][0]) * delta);
                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][2] * B_const[3][0] * (12 * pow(delta, 3) + 3 * (sub_const[1][2] + sub_const[3][0]) * pow(delta, 2) + (sub_const[1][2] * sub_const[3][0]) * delta);

                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][0] * B_const[3][0] * (-12 * pow(delta, 3) + 3 * (sub_const[2][0] - sub_const[3][0]) * pow(delta, 2) + (sub_const[2][0] * sub_const[3][0]) * delta);
                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][1] * B_const[3][0] * (-12 * pow(delta, 3) + 3 * (sub_const[2][1] - sub_const[3][0]) * pow(delta, 2) + (sub_const[2][1] * sub_const[3][0]) * delta);
                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][2] * B_const[3][0] * (-12 * pow(delta, 3) + 3 * (sub_const[2][2] - sub_const[3][0]) * pow(delta, 2) + (sub_const[2][2] * sub_const[3][0]) * delta);
            }
        }
        for (int dim = 0; dim < 0; dim++)
        {
            for (int k = 3; k < u_knot_temp.size() - 4; k++)
            {
                auto u = u_knot_temp;
                vector<vector<double>> B_const;
                B_const.resize(4);
                B_const[0].resize(1);
                B_const[1].resize(3);
                B_const[2].resize(3);
                B_const[3].resize(1);
                B_const[0][0] = 1.0 / ((u[k + 1] - u[k]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k - 2]));

                B_const[1][0] = 1.0 / ((u[k + 1] - u[k - 2]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[1][1] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[1][2] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[2][0] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 1] - u[k - 1]) * (u[k + 1] - u[k]));

                B_const[2][1] = 1.0 / ((u[k + 2] - u[k - 1]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[2][2] = 1.0 / ((u[k + 3] - u[k]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));

                B_const[3][0] = 1.0 / ((u[k + 3] - u[k]) * (u[k + 2] - u[k]) * (u[k + 1] - u[k]));



                double delta = u[k + 1] - u[k];
                Vec3 p0 = P[k - 3][dim];
                Vec3 p1 = P[k - 2][dim];
                Vec3 p2 = P[k - 1][dim];
                Vec3 p3 = P[k][dim];
                s += (p0[0] * p0[0] + p0[1] * p0[1] + p0[2] * p0[2]) * B_const[0][0] * B_const[0][0] * 36 * delta;

                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][0] * 36 * delta;
                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][1] * B_const[1][1] * 36 * delta;
                s += (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][2] * B_const[1][2] * 36 * delta;
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][1] * 36 * delta;
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][0] * B_const[1][2] * 36 * delta;
                s += 2 * (p1[0] * p1[0] + p1[1] * p1[1] + p1[2] * p1[2]) * B_const[1][1] * B_const[1][2] * 36 * delta;

                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][0] * 36 * delta;
                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][1] * B_const[2][1] * 36 * delta;
                s += (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][2] * B_const[2][2] * 36 * delta;
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][1] * 36 * delta;
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][0] * B_const[2][2] * 36 * delta;
                s += 2 * (p2[0] * p2[0] + p2[1] * p2[1] + p2[2] * p2[2]) * B_const[2][1] * B_const[2][2] * 36 * delta;

                s += (p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2]) * B_const[3][0] * B_const[3][0] * 36 * delta;

                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][0] * (-36) * delta;
                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][1] * (-36) * delta;
                s += 2 * (p0[0] * p1[0] + p0[1] * p1[1] + p0[2] * p1[2]) * B_const[0][0] * B_const[1][2] * (-36) * delta;

                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][0] * 36 * delta;
                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][1] * 36 * delta;
                s += 2 * (p0[0] * p2[0] + p0[1] * p2[1] + p0[2] * p2[2]) * B_const[0][0] * B_const[2][2] * 36 * delta;

                s += 2 * (p0[0] * p3[0] + p0[1] * p3[1] + p0[2] * p3[2]) * B_const[0][0] * B_const[3][0] * (-36) * delta;

                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][0] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][1] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][0] * B_const[2][2] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][0] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][1] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][1] * B_const[2][2] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][0] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][1] * (-36) * delta;
                s += 2 * (p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2]) * B_const[1][2] * B_const[2][2] * (-36) * delta;

                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][0] * B_const[3][0] * 36 * delta;
                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][1] * B_const[3][0] * 36 * delta;
                s += 2 * (p1[0] * p3[0] + p1[1] * p3[1] + p1[2] * p3[2]) * B_const[1][2] * B_const[3][0] * 36 * delta;

                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][0] * B_const[3][0] * (-36) * delta;
                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][1] * B_const[3][0] * (-36) * delta;
                s += 2 * (p2[0] * p3[0] + p2[1] * p3[1] + p2[2] * p3[2]) * B_const[2][2] * B_const[3][0] * (-36) * delta;
            }
        }
        //s *= 1e-4;
    }
    scalar_t s1 = 0;
    /**/
    for (int dim_row = 0; dim_row < curve_order_num; dim_row++)
    {
        Vec3 pi0; pi0 = P[dim_row][0];
        Vec3 pi1; pi1 = P[dim_row][1];

        Vec3 pii0; pii0 = P[dim_row + 1][0];
        Vec3 pii1; pii1 = P[dim_row + 1][1];

        Vec3 ruls = (pii1 - pii0) - (pi1 - pi0);
        s1 += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
    }

    vector<double> res = uni_para(5);
    vector<vector<Vec3>> curve_res;
    curve_res.resize(2);

    for (int i = 0; i < res.size(); i++)
    {
        int k = u_k_position(res[i], u_knot_temp);
        vector<double> u_basis_temp;
        calc_basis_fuction(res[i], k, poly_degree, u_basis_temp);
        Vec3 bspline_surface_value1;
        bspline_surface_value1.setZero();
        Vec3 bspline_surface_value2;
        bspline_surface_value2.setZero();
        if (k == 0)
        {

            bspline_surface_value1 = P[0][0];

            bspline_surface_value2 = P[0][1];
        }
        else if (k == u_knot_temp.size() - 1)
        {
            bspline_surface_value1 = P[spline_v.size() / 6 - 1][0];

            bspline_surface_value2 = P[spline_v.size() / 6 - 1][1];
        }
        else
        {
            for (int j = 0; j < u_basis_temp.size(); j++)
            {
                bspline_surface_value1[0] += u_basis_temp[j] * P[k - 3 + j][0][0];
                bspline_surface_value1[1] += u_basis_temp[j] * P[k - 3 + j][0][1];
                bspline_surface_value1[2] += u_basis_temp[j] * P[k - 3 + j][0][2];
                bspline_surface_value2[0] += u_basis_temp[j] * P[k - 3 + j][1][0];
                bspline_surface_value2[1] += u_basis_temp[j] * P[k - 3 + j][1][1];
                bspline_surface_value2[2] += u_basis_temp[j] * P[k - 3 + j][1][2];
            }
        }
        curve_res[0].push_back(bspline_surface_value1);
        curve_res[1].push_back(bspline_surface_value2);

    }

    for (int dim_row = 0; dim_row < curve_res[0].size() - 3; dim_row++)
    {
        Vec3 pi0 = curve_res[0][dim_row];
        Vec3 pi1 = curve_res[1][dim_row];

        Vec3 pii0; pii0 = curve_res[0][dim_row + 3];
        Vec3 pii1; pii1 = curve_res[1][dim_row + 3];

        Vec3 ruls = (pii1 - pii0) - (pi1 - pi0);
        s1 += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
    }
    for (int dim_col = 0; dim_col < 2; dim_col++)
    {
        for (int dim_row = 0; dim_row < curve_res[0].size() - 4; dim_row++)
        {
            Vec3 pi0 = curve_res[dim_col][dim_row];
            Vec3 pi1 = curve_res[dim_col][dim_row + 2];

            Vec3 pii0; pii0 = curve_res[dim_col][dim_row + 4];
            Vec3 ruls = pi0 - pi1 + pii0 - pi1;
            s += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
        }
    }


    scalar_t s_total = u_smooth_weight * s + v_smooth_weight * s1;
    smooth_value = s_total.value().value();
    cout << "smooth_value:=" << smooth_value << endl;
    return smooth_value;
}

double BSplineFitting::line_search(VectorXd& BSpline, Eigen::VectorXd& grad, Eigen::VectorXd& decrease_direction, double c_parm, double beta, double data_term)
{
    double step_size = 1.0;// Initial step size
    int iter_num = 0;
    while (iter_num < 10)
    {
        Eigen::VectorXd new_spline = BSpline + step_size * decrease_direction;
        double data_term_new = calc_data_term_energy(new_spline);
        std::cout << endl;
        std::cout << "setp: " << step_size << endl;
        std::cout << "new_total: " << data_term_new << ", data_value_new: " << data_term_new << ", smooth_value_new: " << 0 << endl;
        std::cout << "old_total: " << data_term << ", data_value: " << data_term << ", smooth_value: " << 0 << endl;

        if (data_term_new > data_term + c_parm * step_size * decrease_direction.transpose() * grad)
        {
            step_size *= beta;
        }
        else
        {
            return step_size;
            break;
        }

        iter_num++;
        if (iter_num == 10)
        {
            step_size = 0;
            return step_size;
        }
    }
}

void BSplineFitting::optimization(Mesh& mesh, Eigen::VectorXd& BSpline)
{
    //设置曲线的节点向量
    u_knot_temp.clear();
    u_knot_temp.shrink_to_fit();
    int KnotNum = BSpline.size() / 6 + poly_degree + 1;
    u_knot_temp.resize(KnotNum);

    int no_repeat_interval = BSpline.size() / 6 - poly_degree;//原本错误
    for (int i = 0; i < KnotNum; i++)
    {
        if (i < poly_degree + 1)
        {
            u_knot_temp[i] = 0;
        }
        else if (i > KnotNum - poly_degree - 1)
        {
            u_knot_temp[i] = no_repeat_interval;
        }
        else
        {
            u_knot_temp[i] = i - poly_degree;
        }

        //cout << u_knot_temp[i] << " ";
    }

    //初始化 细分的控制点及节点向量
    divide_ctrpoint.resize(10);//最大十次细分
    divide_u_knot.resize(10);
    divide_basis.resize(10);
    for (int i = 0; i < 10; i++)
    {
        divide_ctrpoint[i].resize(pow(2, i));
        divide_u_knot[i].resize(pow(2, i));
        divide_basis[i].resize(pow(2, i));
    }
    divide_ctrpoint[0][0] = BSpline;
    divide_u_knot[0][0] = u_knot_temp;
    divide_basis[0][0].resize(BSpline.size() / 6, BSpline.size() / 6);
    divide_basis[0][0].setIdentity();

    for (int i = 0; i < divide_basis[0][0].rows(); i++)
    {
        for (int j = 0; j < divide_basis[0][0].cols(); j++)
        {
            cout << divide_basis[0][0].coeffRef(i, j) << " ";
        }
        cout << endl;
    }

    //test
    /*
    u_knot_temp.resize(11);
    u_knot_temp[0] = 0;
    u_knot_temp[1] = 0;
    u_knot_temp[2] = 0;
    u_knot_temp[3] = 0;

    u_knot_temp[4] = 0.25;
    u_knot_temp[5] = 0.5;
    u_knot_temp[6] = 0.75;

    u_knot_temp[7] = 1;
    u_knot_temp[8] = 1;
    u_knot_temp[9] = 1;
    u_knot_temp[10] = 1;
    */
    divide();
    divide();
    //divide();
    //
    /*
    Eigen::MatrixXd ctr1(BSpline.size() / 6, 3);
    for (int i = 0; i < ctr1.rows(); i++)
    {
        for (int j = 0; j < ctr1.cols(); j++)
        {
            ctr1(i, j) = BSpline[3 * i + j];

        }
        cout << endl;
    }

    Eigen::MatrixXd ctr2(divide_ctrpoint[2][0].size() / 6, 3);
    for (int i = 0; i < ctr2.rows(); i++)
    {
        for (int j = 0; j < ctr2.cols(); j++)
        {
            ctr2(i, j) = divide_ctrpoint[2][0][3 * i + j];
            cout << ctr2.coeffRef(i, j) << " ";

        }
        cout << endl;
    }


    Eigen::MatrixXd test1 = divide_basis[2][0] * ctr1;

    auto dif = test1 - ctr2;
    cout << endl << endl;
    for (int i = 0; i < test1.rows(); i++)
    {
        for (int j = 0; j < test1.cols(); j++)
        {
            cout << test1.coeffRef(i, j) << " ";

        }
        cout << endl;
    }*/

    interval_up_bound = u_knot_temp.back();

    control_number = BSpline.size() / 6;


    BSpline_initialization();
    double stop_threshold = 1e-6;
    int iter_num = 0;

    //BSpline_surface_viewer_2(BSpline, 50, 0);

    Eigen::VectorXd new_BSpline;
    while (true)
    {

        Eigen::VectorXd data_grad;
        Eigen::MatrixXd data_hassien;
        Eigen::VectorXd smooth_grad;
        Eigen::MatrixXd smooth_hassien;
        Eigen::VectorXd decrease_direction;
        Eigen::VectorXd total_grad;
        Eigen::MatrixXd total_hassien;
        scalar_t data_value; scalar_t smooth_value;

        /*
        calc_objective_grad_hessien(BSpline, data_value, data_grad, data_hassien, decrease_direction, smooth_value, smooth_grad, smooth_hassien, total_grad, total_hassien);
        double data_term = data_value.value().value();
        double step = line_search(BSpline, total_grad, decrease_direction, 1e-4, 0.5, data_term);

        new_BSpline = BSpline + step * decrease_direction;
        double termination_condition = (new_BSpline - BSpline).squaredNorm();
        BSpline = new_BSpline;

        iter_num++;


        BSpline_surface_viewer_2(BSpline, 50, iter_num);

        if (termination_condition < stop_threshold || iter_num > 100)
        {
            break;
        }*/
    }
}

void BSplineFitting::test()
{

    /*
    double u_param = 0.5;
    std::vector<double> basis_func;
    u_knot_temp = { 0,0,0,0,1,2,3,4,5,6,6,6,6 };
    auto it = std::upper_bound(u_knot_temp.begin(), u_knot_temp.end(), u_param);
    int interval_order = std::distance(u_knot_temp.begin(), it) - 1;
    calc_basis_fuction(u_param, interval_order, 3, basis_func);

    Handle(Geom_BSplineCurve) curve = nullptr;
    int ctr_num = 9;
    vector<CtrInfo> ctr_point;
    for (int k = 0; k < ctr_num; k++)
    {
        CtrInfo ci;
        ci.position = { Standard_Real(k + 0.5),Standard_Real(k),Standard_Real(2 * k) };
        ci.weight = 1;
        ctr_point.push_back(ci);
    }
    create_BSpline_curve(ctr_point, curve);

    gp_Pnt p0;
    curve->D0(u_param, p0);
    cout << p0.X()<<" " << p0.Y() << " " << p0.Z() << " " << endl;

    Standard_Real px = 0;
    Standard_Real py = 0;
    Standard_Real pz = 0;
    for (int i = 0; i < basis_func.size(); i++)
    {
        px += basis_func[i] * ctr_point[interval_order - poly_degree + i].position.X();
        py += basis_func[i] * ctr_point[interval_order - poly_degree + i].position.Y();
        pz += basis_func[i] * ctr_point[interval_order - poly_degree + i].position.Z();
    }
    cout << px << " " << py << " " << pz << " " << endl;*/
}

void BSplineFitting::divide()
{
    cout << "第" << divide_num + 1 << "次细分：" << endl;
    for (int l = 0; l < pow(2, divide_num); l++)
    {
        double u = (divide_u_knot[divide_num][l][int(divide_u_knot[divide_num][l].size() / 2)] + divide_u_knot[divide_num][l][int(divide_u_knot[divide_num][l].size() / 2) - 1]) / 2;
        auto BSpline = divide_ctrpoint[divide_num][l];
        auto u_knot_temp1 = divide_u_knot[divide_num][l];
        auto temp_ctrpoint = BSpline;
        int k = u_k_position(u, u_knot_temp1);
        cout << "u=" << u << ";k=" << k << endl;
        vector<vector<double>> a;
        a.resize(poly_degree);
        for (int i = 0; i < poly_degree; i++)
        {
            a[i].resize(poly_degree - i);
            for (int j = 0; j < a[i].size(); j++)
            {
                a[i][j] = (u - u_knot_temp1[j + k - poly_degree + i + 1]) / (u_knot_temp1[j + k + 1] - u_knot_temp1[j + k - poly_degree + i + 1]);
                //cout << a[i][j] << " ";
            }
            //cout << endl;
        }

        vector<vector<vector<double>>> ctrpoint_new;
        ctrpoint_new.resize(poly_degree + 1);
        for (int i = 0; i < poly_degree + 1; i++)
        {
            ctrpoint_new[i].resize(poly_degree - i + 1);
            if (i == 0)
            {
                for (int j = 0; j < ctrpoint_new[i].size(); j++)
                {
                    ctrpoint_new[i][j].resize(3);
                    ctrpoint_new[i][j][0] = BSpline[3 * (k - poly_degree + j) + 0];
                    ctrpoint_new[i][j][1] = BSpline[3 * (k - poly_degree + j) + 1];
                    ctrpoint_new[i][j][2] = BSpline[3 * (k - poly_degree + j) + 2];
                    //cout << "P_" << j + k - poly_degree << "," << i << "=" << endl;
                }
            }
            else
            {
                for (int j = 0; j < ctrpoint_new[i].size(); j++)
                {
                    ctrpoint_new[i][j].resize(3);
                    ctrpoint_new[i][j][0] = (1 - a[i - 1][j]) * ctrpoint_new[i - 1][j][0] + a[i - 1][j] * ctrpoint_new[i - 1][j + 1][0];
                    ctrpoint_new[i][j][1] = (1 - a[i - 1][j]) * ctrpoint_new[i - 1][j][1] + a[i - 1][j] * ctrpoint_new[i - 1][j + 1][1];
                    ctrpoint_new[i][j][2] = (1 - a[i - 1][j]) * ctrpoint_new[i - 1][j][2] + a[i - 1][j] * ctrpoint_new[i - 1][j + 1][2];
                    //cout << "P_" << j + k - poly_degree << "," << i << "=" << (1 - a[i - 1][j])<<"*"<<"P_"<< j + k - poly_degree << "," << i-1 <<"+" << (a[i - 1][j]) << "*" << "P_" << j+1 + k - poly_degree << "," << i - 1 << endl;
                }
            }
        }
        vector<vector<vector<double>>> ctrpoint_new2;
        ctrpoint_new2.resize(poly_degree + 1);
        for (int i = 0; i < poly_degree + 1; i++)
        {
            ctrpoint_new2[i].resize(poly_degree - i + 1);
            if (i == 0)
            {
                for (int j = 0; j < ctrpoint_new2[i].size(); j++)
                {
                    ctrpoint_new2[i][j].resize(3);
                    ctrpoint_new2[i][j][0] = BSpline[3 * (k - poly_degree + j) + BSpline.size() / 2 + 0];
                    ctrpoint_new2[i][j][1] = BSpline[3 * (k - poly_degree + j) + BSpline.size() / 2 + 1];
                    ctrpoint_new2[i][j][2] = BSpline[3 * (k - poly_degree + j) + BSpline.size() / 2 + 2];
                    //cout << "P_" << j + k - poly_degree << "," << i << "=" << endl;
                }
            }
            else
            {
                for (int j = 0; j < ctrpoint_new2[i].size(); j++)
                {
                    ctrpoint_new2[i][j].resize(3);
                    ctrpoint_new2[i][j][0] = (1 - a[i - 1][j]) * ctrpoint_new2[i - 1][j][0] + a[i - 1][j] * ctrpoint_new2[i - 1][j + 1][0];
                    ctrpoint_new2[i][j][1] = (1 - a[i - 1][j]) * ctrpoint_new2[i - 1][j][1] + a[i - 1][j] * ctrpoint_new2[i - 1][j + 1][1];
                    ctrpoint_new2[i][j][2] = (1 - a[i - 1][j]) * ctrpoint_new2[i - 1][j][2] + a[i - 1][j] * ctrpoint_new2[i - 1][j + 1][2];
                    //cout << "P_" << j + k - poly_degree << "," << i << "=" << (1 - a[i - 1][j])<<"*"<<"P_"<< j + k - poly_degree << "," << i-1 <<"+" << (a[i - 1][j]) << "*" << "P_" << j+1 + k - poly_degree << "," << i - 1 << endl;
                }
            }
        }

        temp_ctrpoint.resize(BSpline.size() + 3 * 2 * (poly_degree + 1));

        divide_ctrpoint[divide_num + 1][2 * l].resize(3 * 2 * (k + 1));
        divide_ctrpoint[divide_num + 1][2 * l + 1].resize(3 * 2 * (BSpline.size() / 6 - k + poly_degree));
        Eigen::MatrixXd A1(k + 1, BSpline.size() / 6);
        Eigen::MatrixXd A2(BSpline.size() / 6 - k + poly_degree, BSpline.size() / 6);
        A1.setZero();
        A2.setZero();

        for (int i = 0; i < k - poly_degree + 1; i++)
        {
            A1(i, i) = 1;
        }

        A1(k - poly_degree + 1, k - poly_degree) = 1 - a[0][0];
        A1(k - poly_degree + 1, k - poly_degree + 1) = a[0][0];
        A1(k - poly_degree + 2, k - poly_degree) = (1 - a[0][0]) * (1 - a[1][0]);
        A1(k - poly_degree + 2, k - poly_degree + 1) = (a[0][0]) * (1 - a[1][0]) + (1 - a[0][1]) * a[1][0];
        A1(k - poly_degree + 2, k - poly_degree + 2) = (a[0][1]) * a[1][0];
        A1(k - poly_degree + 3, k - poly_degree) = (1 - a[0][0]) * (1 - a[1][0]) * (1 - a[2][0]);
        A1(k - poly_degree + 3, k - poly_degree + 1) = (a[0][0]) * (1 - a[1][0]) * (1 - a[2][0]) + (1 - a[0][1]) * a[1][0] * (1 - a[2][0]) + (1 - a[0][1]) * (1 - a[1][1]) * a[2][0];
        A1(k - poly_degree + 3, k - poly_degree + 2) = (a[0][1]) * a[1][0] * (1 - a[2][0]) + (a[0][1]) * (1 - a[1][1]) * (a[2][0]) + (1 - a[0][2]) * (a[1][1]) * a[2][0];
        A1(k - poly_degree + 3, k - poly_degree + 3) = (a[0][2]) * (a[1][1]) * (a[2][0]);



        divide_basis[divide_num + 1][2 * l] = A1 * divide_basis[divide_num][l];

        for (int i = BSpline.size() / 6 - k + poly_degree - 1; i >= poly_degree; i--)
        {
            A2(i, i + k - poly_degree) = 1;
        }
        A2(2, k) = a[0][2];
        A2(2, k - 1) = 1 - a[0][2];
        A2(1, k) = a[0][2] * a[1][1];
        A2(1, k - 1) = (1 - a[0][2]) * a[1][1] + (a[0][1]) * (1 - a[1][1]);
        A2(1, k - 2) = (1 - a[0][1]) * (1 - a[1][1]);
        A2(0, k) = a[0][2] * a[1][1] * a[2][0];
        A2(0, k - 1) = A2(1, k - 1) * a[2][0] + a[0][1] * a[1][0] * (1 - a[2][0]);
        A2(0, k - 2) = (1 - a[0][1]) * (1 - a[1][1]) * a[2][0] + (1 - a[0][1]) * a[1][0] * (1 - a[2][0]) + a[0][0] * (1 - a[1][0]) * (1 - a[2][0]);
        A2(0, k - 3) = (1 - a[0][0]) * (1 - a[1][0]) * (1 - a[2][0]);

        divide_basis[divide_num + 1][2 * l + 1] = A2 * divide_basis[divide_num][l];
        for (int i = 0; i < (k - poly_degree); i++)
        {
            temp_ctrpoint[3 * i] = BSpline[3 * i];
            temp_ctrpoint[3 * i + 1] = BSpline[3 * i + 1];
            temp_ctrpoint[3 * i + 2] = BSpline[3 * i + 2];

            divide_ctrpoint[divide_num + 1][2 * l][3 * i] = BSpline[3 * i];
            divide_ctrpoint[divide_num + 1][2 * l][3 * i + 1] = BSpline[3 * i + 1];
            divide_ctrpoint[divide_num + 1][2 * l][3 * i + 2] = BSpline[3 * i + 2];
            //cout << i << endl;
        }
        //cout << endl;
        for (int i = 0; i < (poly_degree + 1); i++)
        {
            divide_ctrpoint[divide_num + 1][2 * l][3 * (k - poly_degree + i)] = ctrpoint_new[i][0][0];
            divide_ctrpoint[divide_num + 1][2 * l][3 * (k - poly_degree + i) + 1] = ctrpoint_new[i][0][1];
            divide_ctrpoint[divide_num + 1][2 * l][3 * (k - poly_degree + i) + 2] = ctrpoint_new[i][0][2];
            temp_ctrpoint[3 * (k - poly_degree + i)] = ctrpoint_new[i][0][0];
            temp_ctrpoint[3 * (k - poly_degree + i) + 1] = ctrpoint_new[i][0][1];
            temp_ctrpoint[3 * (k - poly_degree + i) + 2] = ctrpoint_new[i][0][2];
            //cout << k - poly_degree + i << endl;
        }
        //cout << endl;
        for (int i = poly_degree; i >= 0; i--)
        {
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (poly_degree - i)] = ctrpoint_new[i][ctrpoint_new[i].size() - 1][0];
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (poly_degree - i) + 1] = ctrpoint_new[i][ctrpoint_new[i].size() - 1][1];
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (poly_degree - i) + 2] = ctrpoint_new[i][ctrpoint_new[i].size() - 1][2];
            //cout << ( poly_degree - i) << endl;
        }
        //cout << endl;
        for (int i = k + poly_degree + 2; i < temp_ctrpoint.size() / 6; i++)
        {
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (i - k - 1)] = BSpline[3 * (i - (poly_degree + 1))];
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (i - k - 1) + 1] = BSpline[3 * (i - (poly_degree + 1)) + 1];
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (i - k - 1) + 2] = BSpline[3 * (i - (poly_degree + 1)) + 2];
            //cout << (i - k - 1) << endl;
        }

        //cout << endl;
        for (int i = 0; i < (k - poly_degree); i++)
        {
            divide_ctrpoint[divide_num + 1][2 * l][3 * i + divide_ctrpoint[divide_num + 1][2 * l].size() / 2] = BSpline[3 * i + BSpline.size() / 2];
            divide_ctrpoint[divide_num + 1][2 * l][3 * i + 1 + divide_ctrpoint[divide_num + 1][2 * l].size() / 2] = BSpline[3 * i + 1 + BSpline.size() / 2];
            divide_ctrpoint[divide_num + 1][2 * l][3 * i + 2 + divide_ctrpoint[divide_num + 1][2 * l].size() / 2] = BSpline[3 * i + 2 + BSpline.size() / 2];
            //cout << i + divide_ctrpoint[divide_num + 1][2 * l].size() / 6 << endl;
        }
        //cout << endl;
        for (int i = 0; i < (poly_degree + 1); i++)
        {
            divide_ctrpoint[divide_num + 1][2 * l][3 * (k - poly_degree + i) + divide_ctrpoint[divide_num + 1][2 * l].size() / 2] = ctrpoint_new2[i][0][0];
            divide_ctrpoint[divide_num + 1][2 * l][3 * (k - poly_degree + i) + 1 + divide_ctrpoint[divide_num + 1][2 * l].size() / 2] = ctrpoint_new2[i][0][1];
            divide_ctrpoint[divide_num + 1][2 * l][3 * (k - poly_degree + i) + 2 + divide_ctrpoint[divide_num + 1][2 * l].size() / 2] = ctrpoint_new2[i][0][2];
            //cout << k - poly_degree + i + divide_ctrpoint[divide_num + 1][2 * l].size() / 6 << endl;
        }
        //cout << endl;
        for (int i = poly_degree; i >= 0; i--)
        {
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (poly_degree - i) + divide_ctrpoint[divide_num + 1][2 * l + 1].size() / 2] = ctrpoint_new2[i][ctrpoint_new2[i].size() - 1][0];
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (poly_degree - i) + 1 + divide_ctrpoint[divide_num + 1][2 * l + 1].size() / 2] = ctrpoint_new2[i][ctrpoint_new2[i].size() - 1][1];
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (poly_degree - i) + 2 + divide_ctrpoint[divide_num + 1][2 * l + 1].size() / 2] = ctrpoint_new2[i][ctrpoint_new2[i].size() - 1][2];
            //cout << ( poly_degree - i) + divide_ctrpoint[divide_num + 1][2 * l + 1].size() / 6 << endl;
        }
        //cout << endl;
        for (int i = k + poly_degree + 2; i < temp_ctrpoint.size() / 6; i++)
        {
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (i - k - 1) + divide_ctrpoint[divide_num + 1][2 * l + 1].size() / 2] = BSpline[3 * (i - (poly_degree + 1)) + BSpline.size() / 2];
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (i - k - 1) + 1 + divide_ctrpoint[divide_num + 1][2 * l + 1].size() / 2] = BSpline[3 * (i - (poly_degree + 1)) + 1 + BSpline.size() / 2];
            divide_ctrpoint[divide_num + 1][2 * l + 1][3 * (i - k - 1) + 2 + divide_ctrpoint[divide_num + 1][2 * l + 1].size() / 2] = BSpline[3 * (i - (poly_degree + 1)) + 2 + BSpline.size() / 2];
            //cout << (i-k-1) + divide_ctrpoint[divide_num + 1][2 * l + 1].size() / 6 << endl;
        }

        //BSpline = temp_ctrpoint;
        auto knot_temp = u_knot_temp1;
        divide_u_knot[divide_num + 1][2 * l].resize(k + poly_degree + 2);
        divide_u_knot[divide_num + 1][2 * l + 1].resize(BSpline.size() / 6 - k + 2 * (poly_degree + 1) - 1);
        knot_temp.resize(u_knot_temp1.size() + 2 * (poly_degree + 1));
        for (int i = 0; i <= k; i++)
        {
            knot_temp[i] = u_knot_temp1[i];
            //cout << knot_temp[i] << " ";
        }
        for (int i = 0; i < 2 * (poly_degree + 1); i++)
        {
            knot_temp[i + k + 1] = u;
            //cout << knot_temp[i + k + 1] << " ";
        }
        for (int i = k + 2 * (poly_degree + 1) + 1; i < knot_temp.size(); i++)
        {
            knot_temp[i] = u_knot_temp1[i - (2 * poly_degree + 2)];
            //cout << knot_temp[i] << " ";
        }
        for (int i = 0; i < divide_u_knot[divide_num + 1][2 * l].size(); i++)
        {
            divide_u_knot[divide_num + 1][2 * l][i] = knot_temp[i];
            //cout << divide_u_knot[divide_num + 1][2 * l][i] << "  ";
        }
        //cout << endl;
        for (int i = 0; i < divide_u_knot[divide_num + 1][2 * l + 1].size(); i++)
        {
            divide_u_knot[divide_num + 1][2 * l + 1][i] = knot_temp[i + divide_u_knot[divide_num + 1][2 * l].size()];
            //cout << divide_u_knot[divide_num + 1][2 * l+1][i] << "  ";
        }
        //cout << endl;
    }
    divide_num++;
}

void BSplineFitting::add_ctr_point(Eigen::VectorXd& BSpline, double u)
{
    int k = u_k_position(u, u_knot_temp);
    vector<double> a(poly_degree, 0);
    for (int i = 0; i < poly_degree; i++)
    {
        a[i] = (u - u_knot_temp[i + k - poly_degree + 1]) / (u_knot_temp[i + k + 1] - u_knot_temp[i + k - poly_degree + 1]);
    }

    auto ctr_temp = BSpline;
    BSpline.resize(BSpline.size() + 6);
    for (int i = 0; i <= k - poly_degree; i++)
    {
        BSpline[3 * i] = ctr_temp[3 * i];
        BSpline[3 * i + 1] = ctr_temp[3 * i + 1];
        BSpline[3 * i + 2] = ctr_temp[3 * i + 2];

        BSpline[3 * i + BSpline.size() / 2] = ctr_temp[3 * i + ctr_temp.size() / 2];
        BSpline[3 * i + 1 + BSpline.size() / 2] = ctr_temp[3 * i + 1 + ctr_temp.size() / 2];
        BSpline[3 * i + 2 + BSpline.size() / 2] = ctr_temp[3 * i + 2 + ctr_temp.size() / 2];
    }

    for (int i = k - poly_degree + 1; i <= k; i++)
    {
        BSpline[3 * i] = (1 - a[i - (k - poly_degree + 1)]) * ctr_temp[3 * (i - 1)] + a[i - (k - poly_degree + 1)] * ctr_temp[3 * i];
        BSpline[3 * i + 1] = (1 - a[i - (k - poly_degree + 1)]) * ctr_temp[3 * (i - 1) + 1] + a[i - (k - poly_degree + 1)] * ctr_temp[3 * i + 1];
        BSpline[3 * i + 2] = (1 - a[i - (k - poly_degree + 1)]) * ctr_temp[3 * (i - 1) + 2] + a[i - (k - poly_degree + 1)] * ctr_temp[3 * i + 2];

        BSpline[3 * i + BSpline.size() / 2] = (1 - a[i - (k - poly_degree + 1)]) * ctr_temp[3 * (i - 1) + ctr_temp.size() / 2] + a[i - (k - poly_degree + 1)] * ctr_temp[3 * i + ctr_temp.size() / 2];
        BSpline[3 * i + 1 + BSpline.size() / 2] = (1 - a[i - (k - poly_degree + 1)]) * ctr_temp[3 * (i - 1) + 1 + ctr_temp.size() / 2] + a[i - (k - poly_degree + 1)] * ctr_temp[3 * i + 1 + ctr_temp.size() / 2];
        BSpline[3 * i + 2 + BSpline.size() / 2] = (1 - a[i - (k - poly_degree + 1)]) * ctr_temp[3 * (i - 1) + 2 + ctr_temp.size() / 2] + a[i - (k - poly_degree + 1)] * ctr_temp[3 * i + 2 + ctr_temp.size() / 2];
    }

    for (int i = k + 1; i < BSpline.size() / 6; i++)
    {
        BSpline[3 * i] = ctr_temp[3 * (i - 1)];
        BSpline[3 * i + 1] = ctr_temp[3 * (i - 1) + 1];
        BSpline[3 * i + 2] = ctr_temp[3 * (i - 1) + 2];

        BSpline[3 * i + BSpline.size() / 2] = ctr_temp[3 * (i - 1) + ctr_temp.size() / 2];
        BSpline[3 * i + 1 + BSpline.size() / 2] = ctr_temp[3 * (i - 1) + 1 + ctr_temp.size() / 2];
        BSpline[3 * i + 2 + BSpline.size() / 2] = ctr_temp[3 * (i - 1) + 2 + ctr_temp.size() / 2];
    }

    auto knot_temp = u_knot_temp;
    u_knot_temp.clear();
    u_knot_temp.shrink_to_fit();
    u_knot_temp.resize(knot_temp.size() + 1);
    for (int i = 0; i <= k; i++)
    {
        u_knot_temp[i] = knot_temp[i];
    }
    u_knot_temp[k + 1] = u;
    for (int i = k + 2; i < u_knot_temp.size(); i++)
    {
        u_knot_temp[i] = knot_temp[i - 1];
    }
}



void BSplineFitting::CCD_initialization(Data& spline)
{


    VectorXd bspline;
    bspline.resize(spline.rows() * 3);
    for (int i = 0; i < spline.rows(); i++)
    {
        bspline[3 * i] = spline.coeffRef(i, 0);
        bspline[3 * i + 1] = spline.coeffRef(i, 1);
        bspline[3 * i + 2] = spline.coeffRef(i, 2);
    }


    for (int i = 0; i < mesh.n_faces(); i++)
    {
        cover_area_temp.push_back(i);
    }

    u_knot_temp.clear();
    u_knot_temp.shrink_to_fit();
    int KnotNum = bspline.size() / 6 + poly_degree + 1;
    u_knot_temp.resize(KnotNum);
    int no_repeat_interval = bspline.size() / 6 - poly_degree;//原本错误
    for (int i = 0; i < KnotNum; i++)
    {
        if (i < poly_degree + 1)
        {
            u_knot_temp[i] = 0;
        }
        else if (i > KnotNum - poly_degree - 1)
        {
            u_knot_temp[i] = no_repeat_interval;
        }
        else
        {
            u_knot_temp[i] = i - poly_degree;
        }
    }

    interval_up_bound = u_knot_temp.back();
    control_number = bspline.size() / 6;

    BSpline_initialization();
}


std::vector<double> BSplineFitting::uni_para(int samp_num)
{
    vector<double> res;
    res.resize((u_knot_temp.size() - 3 * 2 - 1) * samp_num + 1);
    for (int i = 0; i < (u_knot_temp.size() - poly_degree * 2 - 1); i++)
    {
        for (int j = 0; j < samp_num; j++)
        {
            res[i * samp_num + j] = (u_knot_temp[i + 1 + poly_degree] * (j)+u_knot_temp[i + poly_degree] * (samp_num - j)) / samp_num;
        }
    }
    res[(u_knot_temp.size() - 3 * 2 - 1) * samp_num] = u_knot_temp.back();
    return res;
}
