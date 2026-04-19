/// 文件名:     LinearPrism.C
/**
  * @ Attention 1: 三棱柱的Jacobian矩阵每个积分点都不一样，所以不能像四面体一样取
  * */



#include "LinearPrism.h"

LinearPrism::LinearPrism(const string& name) : BaseElement<NDIM>(name) {}

/*************************************************************************
 * 析构函数.
 ************************************************************************/
LinearPrism::~LinearPrism() {}

void LinearPrism::buildRecoveryMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id) {
  // TODO: 由用户实现恢复型算法
}

void LinearPrism::buildStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, double T_val) {
  // TODO: 由用户实现刚度矩阵计算
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearPrism");

  /// 取出三棱柱形函数对象
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearPrism");

  /// 取出材料管理器
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  tbox::Pointer<Material> material
      = material_manager->getMaterial(GET_USER_MAT(entity_id));
  /// 取出单元上自由度数目 (三棱柱这里 n_dof 会动态返回 6)
  int n_dof = shape_func->getNumberOfDof();
  ele_mat->resize(n_dof * NDIM, n_dof * NDIM);
  for (int i = 0; i < n_dof * NDIM; ++i) {
    for (int j = 0; j < n_dof * NDIM; ++j) {
      (*ele_mat)(i, j) = 0.0;
    }
  }
  /// 提取积分点、雅可比与体积等所有底层信息
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  /// 与四面体不同，这个地方只需要它的局部积分点信息
  tbox::Array<hier::DoubleVector<NDIM> > local_quad_pnt =
      integrator->getLocalQuadraturePoints();
  //  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出形函数在各个积分点的梯度值
  /// 在Prism里和单纯型完全不同！
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, local_quad_pnt);
  tbox::Array<tbox::Array<double> > moduli = material->getModuli(T_val);

  double a = moduli[0][0];
  double b = moduli[0][1];
  double c = moduli[3][3];
  for (int l = 0; l < num_quad_pnts; ++l){
    double detJ = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
    // 计算当前积分点的体积微元
    double JxW = detJ * weight[l];
    // 将刚度项累加到单元矩阵中
    for (int i = 0; i < n_dof; ++i) {
      for (int j = 0; j < n_dof; ++j) {
        (*ele_mat)(NDIM* i, NDIM* j) +=
            JxW * (a * bas_grad[l][i][0] * bas_grad[l][j][0] +
            c * bas_grad[l][i][1] * bas_grad[l][j][1] +
            c * bas_grad[l][i][2] * bas_grad[l][j][2]);

        (*ele_mat)(NDIM* i + 1, NDIM * j + 1) +=
            JxW * (c * bas_grad[l][i][0] * bas_grad[l][j][0] +
            a * bas_grad[l][i][1] * bas_grad[l][j][1] +
            c * bas_grad[l][i][2] * bas_grad[l][j][2]);

        (*ele_mat)(NDIM* i + 2, NDIM * j + 2) +=
            JxW * (c * bas_grad[l][i][0] * bas_grad[l][j][0] +
            c * bas_grad[l][i][1] * bas_grad[l][j][1] +
            a * bas_grad[l][i][2] * bas_grad[l][j][2]);

        (*ele_mat)(NDIM* i, NDIM* j + 1) +=
            JxW * b * bas_grad[l][i][0] * bas_grad[l][j][1] +
            JxW * c * bas_grad[l][j][0] * bas_grad[l][i][1];

        (*ele_mat)(NDIM* i, NDIM* j + 2) +=
            JxW * b * bas_grad[l][i][0] * bas_grad[l][j][2] +
            JxW * c * bas_grad[l][j][0] * bas_grad[l][i][2];

        (*ele_mat)(NDIM* i + 1, NDIM * j + 2) +=
            JxW * b * bas_grad[l][i][1] * bas_grad[l][j][2] +
            JxW * c * bas_grad[l][j][1] * bas_grad[l][i][2];

        (*ele_mat)(NDIM* i + 1, NDIM * j) +=
            JxW * b * bas_grad[l][i][1] * bas_grad[l][j][0] +
            JxW * c * bas_grad[l][j][1] * bas_grad[l][i][0];

        (*ele_mat)(NDIM* i + 2, NDIM * j) +=
            JxW * b * bas_grad[l][i][2] * bas_grad[l][j][0] +
            JxW * c * bas_grad[l][j][2] * bas_grad[l][i][0];

        (*ele_mat)(NDIM* i + 2, NDIM * j + 1) +=
            JxW * b * bas_grad[l][i][2] * bas_grad[l][j][1] +
            JxW * c * bas_grad[l][j][2] * bas_grad[l][i][1];
      }
    }


  }
}

void LinearPrism::buildMassElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, double T_val) {
  // TODO: 由用户实现质量矩阵计算
}

void LinearPrism::buildDumpElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, double T_val) {
  // TODO: 由用户实现阻尼矩阵计算
}

void LinearPrism::buildElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val) {
  // TODO: 由用户实现单元总矩阵组装
}

void LinearPrism::buildElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    NewmarkData *d_newmark, int entity_id, tbox::Array<double> T_val,
    tbox::Array<double> Tolder_val) {
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearPrism");

  /// 取出三棱柱形函数对象
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearPrism");

  /// 取出材料管理器
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  /// 计算单元平均温度 (注意：三棱柱有 6 个节点，所以除以 6.0)
  double e_Temperature = 0;
  for(int i = 0; i < 6; i++) {
    e_Temperature += T_val[i] / 6.0;
  }
  tbox::Pointer<Material> material
      = material_manager->getMaterial(GET_USER_MAT(entity_id));
  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_vec->resize(n_dof*NDIM);
  for (int j = 0; j < n_dof*NDIM; ++j) {
    (*ele_vec)[j] = 0.0;
  }
  /// 与四面体不同，这个地方只需要它的局部积分点信息
  tbox::Array<hier::DoubleVector<NDIM> > local_quad_pnt =
      integrator->getLocalQuadraturePoints();
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  tbox::Array<double> weight = integrator->getQuadratureWeights();
  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, local_quad_pnt);
  /// 在Prism里和单纯型完全不同！
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, local_quad_pnt);
  //模量矩阵
  tbox::Array<tbox::Array<double> > moduli = material->getModuli(e_Temperature);
  double a = moduli[0][0];
  double b = moduli[0][1];
  //计算热应变
  //应变矩阵=th_stress[1 1 1 0 0 0]
  double alpha=material->getAlpha(e_Temperature);
  double th_stress=alpha*(e_Temperature-293.15);
  double matrix_BTD[n_dof*NDIM];
  for(int j1=0;j1<n_dof*NDIM;j1++)
    matrix_BTD[j1]=0;
  //模量矩阵*热应力=termal_strain*[1 1 1 0 0 0]T
  double thermal_strain=(a+2*b)*th_stress;
  for (int l = 0; l < num_quad_pnts; ++l) {
    double detJ = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
    double JxW = detJ * weight[l];
    for (int i = 0; i < n_dof; ++i) {
      for (int j = 0; j < NDIM; ++j) {
        matrix_BTD[i*NDIM+j] += thermal_strain*bas_grad[l][i][j]*JxW;
      }
    }
  }
  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l){
      for (int j = 0; j < NDIM; ++j)
        (*ele_vec)[i * NDIM + j] += matrix_BTD[i*NDIM+j];
    }
  }

}

void LinearPrism::buildTh_ElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val) {
  // TODO: 由用户实现瞬态热传导矩阵
  /// Done on 2026-04-01
  /// 取出三棱柱积分器对象
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearPrism");

  /// 取出三棱柱形函数对象
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearPrism");

  /// 取出材料管理器
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  /// 计算单元平均温度 (注意：三棱柱有 6 个节点，所以除以 6.0)
  double e_Temperature = 0;
  for(int i = 0; i < 6; i++) {
    e_Temperature += T_val[i] / 6.0;
  }
  tbox::Pointer<Material> material
      = material_manager->getMaterial(GET_USER_MAT(entity_id));
  /// 取出单元上自由度数目 (三棱柱这里 n_dof 会动态返回 6)
  int n_dof = shape_func->getNumberOfDof();
  ele_mat->resize(n_dof, n_dof);
  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      (*ele_mat)(i, j) = 0.0;
    }
  }
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 与四面体不同，这个地方只需要它的局部积分点信息
  tbox::Array<hier::DoubleVector<NDIM> > local_quad_pnt =
      integrator->getLocalQuadraturePoints();
  //  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出形函数在各个积分点的梯度值
  /// 在Prism里和单纯型完全不同！
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, local_quad_pnt);
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, local_quad_pnt);

  /// 根据平均温度获取当前材料的热导率 K
  double density=material->getDensity(e_Temperature);
  double K=material->getK(e_Temperature);
  double Cp=material->getCp(e_Temperature);

  for (int l = 0; l < num_quad_pnts; ++l){
    double detJ = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
    // 计算当前积分点的体积微元
    double JxW = detJ * weight[l];
    // 将刚度项累加到单元矩阵中
    for (int i = 0; i < n_dof; ++i) {
      for (int j = 0; j < n_dof; ++j) {
        (*ele_mat)(i,j) += JxW*(density*Cp*bas_val[l][i]*bas_val[l][j]+
                                (bas_grad[l][i][0]*bas_grad[l][j][0]+
            bas_grad[l][i][1]*bas_grad[l][j][1]+
            bas_grad[l][i][2]*bas_grad[l][j][2])*K*dt);
      }
    }


  }
}

void LinearPrism::buildStaticTh_ElementMatrix(tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
                                              const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val,
                                              tbox::Array<hier::DoubleVector<NDIM> > U_val) {
  // TODO: 由用户实现稳态热传导矩阵
  /// Done on 2026-04-01
  /// 取出三棱柱积分器对象
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearPrism");

  /// 取出三棱柱形函数对象
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearPrism");

  /// 取出材料管理器
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  /// 计算单元平均温度 (注意：三棱柱有 6 个节点，所以除以 6.0)
  double e_Temperature = 0;
  for(int i = 0; i < 6; i++) {
    e_Temperature += T_val[i] / 6.0;
  }
  tbox::Pointer<Material> material
      = material_manager->getMaterial(GET_USER_MAT(entity_id));
  /// 取出单元上自由度数目 (三棱柱这里 n_dof 会动态返回 6)
  int n_dof = shape_func->getNumberOfDof();
  ele_mat->resize(n_dof, n_dof);
  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      (*ele_mat)(i, j) = 0.0;
    }
  }
  /// 提取积分点、雅可比与体积等所有底层信息
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  /// 与四面体不同，这个地方只需要它的局部积分点信息
  tbox::Array<hier::DoubleVector<NDIM> > local_quad_pnt =
      integrator->getLocalQuadraturePoints();
//  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出形函数在各个积分点的梯度值
  /// 在Prism里和单纯型完全不同！
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, local_quad_pnt);
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, local_quad_pnt);
  bool fluid_term = entity_id>100;
  double K=material->getK(e_Temperature);
  double density = material->getDensity(e_Temperature);
  double Cp = material->getCp(e_Temperature);

  for (int l = 0; l < num_quad_pnts; ++l){
    double detJ = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
    // 计算当前积分点的体积微元
    double JxW = detJ * weight[l];
    double K_eff = K;
    /// 只有流体热要处理以下情形
    double u_k = 0.0, v_k = 0.0, w_k = 0.0;
    double tau_T = 0.0;
    if (fluid_term) {
      for (int m = 0; m < n_dof; ++m) {
        u_k += bas_val[l][m] * U_val[m][0]; // U_val 是你传入的已收敛的流场速度
        v_k += bas_val[l][m] * U_val[m][1];
        w_k += bas_val[l][m] * U_val[m][2];
      }
      double U_norm = sqrt(u_k*u_k + v_k*v_k + w_k*w_k);
      // 计算热学 SUPG 稳定化参数 (仅限流体)

      if (fluid_term && U_norm > 1e-12) {
        // 计算特征尺度 u_grad_N_sum (沿流线方向的网格尺度)
        double u_grad_N_sum = 0.0;
        for (int m = 0; m < n_dof; ++m) {
          u_grad_N_sum += fabs(u_k * bas_grad[l][m][0] + v_k * bas_grad[l][m][1] + w_k * bas_grad[l][m][2]);
        }

        if (u_grad_N_sum > 1e-12) {
          double h_supg = 2.0 * U_norm / u_grad_N_sum;
          double alpha = K / (density * Cp); // 热扩散率
          // 热学 SUPG 公式：对流极限 vs 扩散极限
          tau_T = 1.0 / sqrt(pow(2.0 * U_norm / h_supg, 2.0) + pow(4.0 * alpha / (h_supg * h_supg), 2.0));
          // 加上适当的缩放系数防止过度耗散 (类似流体力学)
          double tune_thermal = 0.2;
          tau_T = tau_T * tune_thermal;
        }
        double grad_N_squared_sum = 0.0;
        for (int m = 0; m < n_dof; ++m) {
          double dNm_dx = bas_grad[l][m][0];
          double dNm_dy = bas_grad[l][m][1];
          double dNm_dz = bas_grad[l][m][2];
          grad_N_squared_sum += (dNm_dx * dNm_dx + dNm_dy * dNm_dy + dNm_dz * dNm_dz);
        }

        double h_crosswind = 0.0;
        if (grad_N_squared_sum > 1e-15) {
          // 这个公式会自动过滤掉长边，极其精准地锁定边界层法向 1μm 的极薄厚度！
          h_crosswind = 2.0 / sqrt(grad_N_squared_sum);
        } else {
          h_crosswind = pow(detJ, 1.0 / 3.0); // 仅防万一的退化保护
        }

        // 侧风扩散系数，0.05 是非常稳妥的工程经验值
        double tune_crosswind = 0.05;

        // 🚨 注意这里换成了 h_crosswind ！！！
        double K_art = tune_crosswind * density * Cp * U_norm * h_crosswind;

        // 更新有效热导率
        K_eff = K + K_art;
      }

    }
    // 将刚度项累加到单元矩阵中
    for (int i = 0; i < n_dof; ++i) {
      for (int j = 0; j < n_dof; ++j) {
        double diff_ij = (bas_grad[l][i][0]*bas_grad[l][j][0]+
            bas_grad[l][i][1]*bas_grad[l][j][1]+
            bas_grad[l][i][2]*bas_grad[l][j][2])*K_eff;
        double conv_ij = 0.,supg_ij = 0.;
        if (fluid_term){
          double U_dot_gradNj = u_k * bas_grad[l][j][0]
              + v_k * bas_grad[l][j][1]
              + w_k * bas_grad[l][j][2];
          double U_dot_gradNi = u_k * bas_grad[l][i][0]
              + v_k * bas_grad[l][i][1]
              + w_k * bas_grad[l][i][2];
          conv_ij = density * Cp * bas_val[l][i] * U_dot_gradNj;
          supg_ij = tau_T * (U_dot_gradNi) * (density * Cp * U_dot_gradNj);
        }
        (*ele_mat)(i, j) += JxW*(diff_ij + conv_ij + supg_ij);
      }
    }
  }
}

void LinearPrism::buildTh_ElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    int entity_id, tbox::Array<double> T_val, double e_ThermalSource) {
  // TODO: 由用户实现瞬态热传导右端项
  /// Done on 2026-04-03
  /// 取出三棱柱积分器对象
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearPrism");

  /// 取出三棱柱形函数对象
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearPrism");

  /// 取出材料管理器
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  /// 计算单元平均温度 (注意：三棱柱有 6 个节点，所以除以 6.0)
  double e_Temperature = 0;
  for(int i = 0; i < 6; i++) {
    e_Temperature += T_val[i] / 6.0;
  }
  tbox::Pointer<Material> material
      = material_manager->getMaterial(GET_USER_MAT(entity_id));
  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_vec->resize(n_dof);
  for (int j = 0; j < n_dof; ++j) {
    (*ele_vec)[j] = 0.0;
  }
  /// 与四面体不同，这个地方只需要它的局部积分点信息
  tbox::Array<hier::DoubleVector<NDIM> > local_quad_pnt =
      integrator->getLocalQuadraturePoints();
  tbox::Array<double> weight = integrator->getQuadratureWeights();
  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, local_quad_pnt);
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  //先计算单元质量矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_M = new tbox::Matrix<double>();
  ele_mat_M->resize(n_dof,n_dof);
  double density=material->getDensity(e_Temperature);
  double Cp=material->getCp(e_Temperature);
  for (int l = 0; l < num_quad_pnts; ++l) {
    double detJ = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
    double JxW = detJ * weight[l];
    for (int i = 0; i < n_dof; ++i) {
      for (int j = 0; j < n_dof; ++j) {
        (*ele_mat_M)(i,j) += JxW*density*Cp*(bas_val[l][i]*bas_val[l][j]);
      }
    }
  }
  for (int i1 = 0; i1 < n_dof; ++i1) {
    for (int j = 0; j < n_dof; ++j) {
      (*ele_vec)[i1]+=((*ele_mat_M)(i1, j))*T_val[j];
    }
  }

  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l){
      double detJ = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
      // 计算当前积分点的体积微元
      double JxW = detJ * weight[l];
      (*ele_vec)[i] += JxW * e_ThermalSource * bas_val[l][i]* dt;
    }
  }
}

void LinearPrism::buildStaticTh_ElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    int entity_id, tbox::Array<double> T_val, double e_ThermalSource) {
  // TODO: 由用户实现稳态热传导右端项
  /// Done on 2026-04-03
  /// 取出三棱柱积分器对象
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearPrism");

  /// 取出三棱柱形函数对象
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearPrism");

  /// 取出材料管理器
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  /// 计算单元平均温度 (注意：三棱柱有 6 个节点，所以除以 6.0)
  double e_Temperature = 0;
  for(int i = 0; i < 6; i++) {
    e_Temperature += T_val[i] / 6.0;
  }
  tbox::Pointer<Material> material
      = material_manager->getMaterial(GET_USER_MAT(entity_id));
  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_vec->resize(n_dof);
  for (int j = 0; j < n_dof; ++j) {
    (*ele_vec)[j] = 0.0;
  }
  /// 与四面体不同，这个地方只需要它的局部积分点信息
  tbox::Array<hier::DoubleVector<NDIM> > local_quad_pnt =
      integrator->getLocalQuadraturePoints();
  tbox::Array<double> weight = integrator->getQuadratureWeights();
  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, local_quad_pnt);
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  //先计算单元质量矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_M = new tbox::Matrix<double>();
  ele_mat_M->resize(n_dof,n_dof);
  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l){
      double detJ = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
      // 计算当前积分点的体积微元
      double JxW = detJ * weight[l];
      (*ele_vec)[i] += JxW * e_ThermalSource * bas_val[l][i];
    }
  }

}

void LinearPrism::buildE_ElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val) {
  // TODO: 由用户实现电场矩阵计算
}

void LinearPrism::buildE_ElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    int entity_id, tbox::Array<double> T_val) {
  // TODO: 由用户实现电场右端项计算
}

void LinearPrism::buildInitFluidElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val, double *mu) {
  // TODO: 由用户实现斯托克斯流冷启动矩阵
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
        IntegratorManager<NDIM>::getManager();
    tbox::Pointer<BaseIntegrator<NDIM> > integrator =
        integrator_manager->getIntegrator("LinearPrism");

    tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
        ShapeFunctionManager<NDIM>::getManager();
    tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
        shape_manager->getShapeFunction("LinearPrism");

    /// 2. 取出材料并计算平均温度
    tbox::Pointer<MaterialManager<NDIM> > material_manager =
        MaterialManager<NDIM>::getManager();
    double e_Temperature = 0;
    for(int i = 0; i < 6; i++) {
      e_Temperature += T_val[i] / 6.0;
    }
    tbox::Pointer<Material> material
        = material_manager->getMaterial(GET_USER_MAT(entity_id));
    /// 三棱柱固定有 6 个节点
    /// 每个节点 4 个自由度 (u, v, w, p)，局部矩阵维度 24x24
    int n_dof = shape_func->getNumberOfDof();
    int n_dof_dis = NDIM+1;

    int n_dof_total = n_dof * n_dof_dis;     // 流体单元总自由度: 4 * (3 + 1) = 16

    ele_mat->resize(n_dof_total, n_dof_total);
    for (int i = 0; i < n_dof_total; ++i) {
      for (int j = 0; j < n_dof_total; ++j) {
        (*ele_mat)(i, j) = 0.0;
      }
    }
    double mu_value = (*mu);
    double rho = material->getDensity(e_Temperature);

    /// 取出积分点数目.
    int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
    /// 取出模板单元的面积.
    double volume = integrator->getElementVolume();
    /// 与四面体不同，这个地方只需要它的局部积分点信息
    tbox::Array<hier::DoubleVector<NDIM> > local_quad_pnt =
        integrator->getLocalQuadraturePoints();
    //  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
    tbox::Array<double> weight = integrator->getQuadratureWeights();

    /// 取出形函数在各个积分点的梯度值
    /// 在Prism里和单纯型完全不同！
    tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
        shape_func->gradient(real_vertex, local_quad_pnt);
    tbox::Array<tbox::Array<double> > bas_val =
        shape_func->value(real_vertex, local_quad_pnt);
    for (int l = 0; l < num_quad_pnts; ++l){
      /// 动态计算雅可比行列式 |J|
      double jac = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
      /// 真实的物理体积微元 JxW
      double JxW = jac * weight[l];
      /// 3. Shakib-Hughes 动态各向异性特征长度与 PSPG/SUPG 稳定化
      double u_k = 0.0, v_k = 0.0, w_k = 0.0;
      double U_norm = sqrt(u_k * u_k + v_k * v_k + w_k * w_k);

      double u_grad_N_sum = 0.0;
      double grad_N_squared_sum = 0.0;
      for (int m = 0; m < n_dof; ++m) {
        double dNm_dx = bas_grad[l][m][0];
        double dNm_dy = bas_grad[l][m][1];
        double dNm_dz = bas_grad[l][m][2];

        u_grad_N_sum += fabs(u_k * dNm_dx + v_k * dNm_dy + w_k * dNm_dz);
        grad_N_squared_sum += (dNm_dx * dNm_dx + dNm_dy * dNm_dy + dNm_dz * dNm_dz);
      }

      double h_supg = 0.0;
      double h_pspg = 0.0;
      /// 对于压力连续性方程，它不具备流动方向性偏好，
      /// 使用基于形函数几何梯度的调和平均，自动捕捉三棱柱的“最短边（厚度）”作为防过度耗散的基准尺度
      if (grad_N_squared_sum > 1e-15) {
        h_pspg = 2.0 / sqrt(grad_N_squared_sum);
      } else {
        h_pspg = pow(volume * jac, 1.0 / 3.0); // 极小概率退化保护
      }

      /// PSPG 稳定化参数
      double tau_pspg = 0.01*(h_pspg * h_pspg) / (12.0 * mu_value);

      /// (如果未来用于 N-S Jacobian，这里计算 tau_supg)
      double tau_supg = 0.0;
      if (U_norm > 1e-12 && u_grad_N_sum > 1e-12) {
        h_supg = 2.0 * U_norm / u_grad_N_sum;
        tau_supg = 1.0 / sqrt(pow(2.0 * U_norm / h_supg, 2.0) + pow(4.0 * mu_value / (rho * h_supg * h_supg), 2.0));
      }
      else {
        // 纯扩散极限下，复用 h_pspg 作为特征尺度，避免截断为 0
        tau_supg = (h_pspg * h_pspg) / (4.0 * mu_value / rho);
      }
      // 尝试 0.5, 0.2, 或者 0.1。值越小，中心峰值速度越高（越接近理论抛物线）
      double tune_supg = 0.5e-1;
      tau_supg *= tune_supg;
      for (int i = 0; i < n_dof; ++i){
        int row_u = i * n_dof_dis + 0;
        int row_v = i * n_dof_dis + 1;
        int row_w = i * n_dof_dis + 2;
        int row_p = i * n_dof_dis + 3;
        for (int j = 0; j < n_dof; ++j){
          int col_u = j * n_dof_dis + 0;
          int col_v = j * n_dof_dis + 1;
          int col_w = j * n_dof_dis + 2;
          int col_p = j * n_dof_dis + 3;
          double diff = JxW * mu_value * (bas_grad[l][i][0] * bas_grad[l][j][0] +
              bas_grad[l][i][1] * bas_grad[l][j][1] +
              bas_grad[l][i][2] * bas_grad[l][j][2]);

          double G_x = -JxW * bas_grad[l][i][0] * bas_val[l][j];
          double G_y = -JxW * bas_grad[l][i][1] * bas_val[l][j];
          double G_z = -JxW * bas_grad[l][i][2] * bas_val[l][j];


          double D_x = -JxW * bas_val[l][i] * bas_grad[l][j][0];
          double D_y = -JxW * bas_val[l][i] * bas_grad[l][j][1];
          double D_z = -JxW * bas_val[l][i] * bas_grad[l][j][2];


          double pspg = -JxW * tau_pspg * (bas_grad[l][i][0] * bas_grad[l][j][0] +
              bas_grad[l][i][1] * bas_grad[l][j][1] +
              bas_grad[l][i][2] * bas_grad[l][j][2]);
          // 1. K_mu (对角线上的速度-速度耦合)
          (*ele_mat)(row_u, col_u) += diff;
          (*ele_mat)(row_v, col_v) += diff;
          (*ele_mat)(row_w, col_w) += diff;

          // 2. G (右上角块：速度-压力耦合)
          (*ele_mat)(row_u, col_p) += G_x;
          (*ele_mat)(row_v, col_p) += G_y;
          (*ele_mat)(row_w, col_p) += G_z;

          // 3. D (左下角块：压力-速度耦合，由于乘了-1，它就是 G 的转置)
          (*ele_mat)(row_p, col_u) += D_x;
          (*ele_mat)(row_p, col_v) += D_y;
          (*ele_mat)(row_p, col_w) += D_z;

          // 4. K_pspg (右下角块：压力-压力耦合，拯救主元为0)
          (*ele_mat)(row_p, col_p) += pspg;
        }
      }

    }

}

void LinearPrism::buildFluidJacobianElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,
    int entity_id, tbox::Array<double> T_val,
    tbox::Array<hier::DoubleVector<NDIM> > U_val) {
  // TODO: 由用户实现N-S牛顿雅可比矩阵
  /// Done on 2026-04-14
  /// 取出积分器与形函数对象
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager
      = IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator
      = integrator_manager->getIntegrator("LinearPrism");
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager
      = ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func
      = shape_manager->getShapeFunction("LinearPrism");
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();
  tbox::Pointer<Material> material
      = material_manager->getMaterial(GET_USER_MAT(entity_id));
  double e_Temperature = 0;
  for(int i = 0; i < 6; i++) {
    e_Temperature += T_val[i] / 6.0;
  }
  double mu_value = 1.002e-3;
  double rho = material->getDensity(e_Temperature);
  TBOX_ASSERT(rho>0);
  int n_nodes = shape_func->getNumberOfDof(); // 6
  int n_dof_total = n_nodes * (NDIM + 1);     // 24
  ele_mat->resize(n_dof_total, n_dof_total);
  for (int i = 0; i < n_dof_total; ++i) {
    for (int j = 0; j < n_dof_total; ++j) {
      (*ele_mat)(i, j) = 0.0;
    }
  }
  /// 提取积分点信息
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  double volume = integrator->getElementVolume();
  /// 与四面体不同，这个地方只需要它的局部积分点信息
  tbox::Array<hier::DoubleVector<NDIM> > local_quad_pnt =
      integrator->getLocalQuadraturePoints();
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, local_quad_pnt);
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, local_quad_pnt);
  for (int l = 0; l < num_quad_pnts; ++l){
    /// 动态计算雅可比行列式 |J|
    double jac = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
    /// 真实的物理体积微元 JxW
    double JxW = jac * weight[l];
    /// 3. Shakib-Hughes 动态各向异性特征长度与 PSPG/SUPG 稳定化
    double u_k = 0.0, v_k = 0.0, w_k = 0.0;
    double du_dx = 0.0, du_dy = 0.0, du_dz = 0.0;
    double dv_dx = 0.0, dv_dy = 0.0, dv_dz = 0.0;
    double dw_dx = 0.0, dw_dy = 0.0, dw_dz = 0.0;

    for (int m = 0; m < n_nodes; ++m) {
      double N_m = bas_val[l][m];
      double dNm_dx = bas_grad[l][m][0];
      double dNm_dy = bas_grad[l][m][1];
      double dNm_dz = bas_grad[l][m][2];

      u_k += N_m * U_val[m][0];
      v_k += N_m * U_val[m][1];
      w_k += N_m * U_val[m][2];

      du_dx += dNm_dx * U_val[m][0]; du_dy += dNm_dy * U_val[m][0]; du_dz += dNm_dz * U_val[m][0];
      dv_dx += dNm_dx * U_val[m][1]; dv_dy += dNm_dy * U_val[m][1]; dv_dz += dNm_dz * U_val[m][1];
      dw_dx += dNm_dx * U_val[m][2]; dw_dy += dNm_dy * U_val[m][2]; dw_dz += dNm_dz * U_val[m][2];

    }
    double U_norm = sqrt(u_k * u_k + v_k * v_k + w_k * w_k);
    /// ====================================================================
    /// 第 2 步：拿着完整的流速 (u_k, v_k, w_k)，开一个新循环计算特征尺度参数
    /// ====================================================================
    double u_grad_N_sum = 0.0;
    double grad_N_squared_sum = 0.0;

    for (int m = 0; m < n_nodes; ++m) {
      double dNm_dx = bas_grad[l][m][0];
      double dNm_dy = bas_grad[l][m][1];
      double dNm_dz = bas_grad[l][m][2];

      // 核心：完整的流速向量 与 当前节点形函数梯度的内积 的绝对值
      u_grad_N_sum += fabs(u_k * dNm_dx + v_k * dNm_dy + w_k * dNm_dz);

      // 顺便把 PSPG 需要的纯几何对角线也累加了
      grad_N_squared_sum += (dNm_dx * dNm_dx + dNm_dy * dNm_dy + dNm_dz * dNm_dz);
    }


    double h_supg = 0.0;
    double h_pspg = 0.0;
    /// 对于压力连续性方程，它不具备流动方向性偏好，
    /// 使用基于形函数几何梯度的调和平均，自动捕捉三棱柱的“最短边（厚度）”作为防过度耗散的基准尺度
    if (grad_N_squared_sum > 1e-15) {
      h_pspg = 2.0 / sqrt(grad_N_squared_sum);
    } else {
      h_pspg = pow(volume * jac, 1.0 / 3.0); // 极小概率退化保护
    }

    /// PSPG 稳定化参数
//    double tau_pspg = (h_pspg * h_pspg) / (12.0 * mu_value);

    /// (用于 N-S Jacobian，这里计算 tau_supg)
    double tau_supg = 0.0;
    if (U_norm > 1e-12 && u_grad_N_sum > 1e-12) {
      h_supg = 2.0 * U_norm / u_grad_N_sum;
      tau_supg = 1.0 / sqrt(pow(2.0 * U_norm / h_supg, 2.0) + pow(4.0 * mu_value / (rho * h_supg * h_supg), 2.0));
    }
    else {
      // 纯扩散极限下，复用 h_pspg 作为特征尺度，避免截断为 0
      tau_supg = (h_pspg * h_pspg) / (4.0 * mu_value / rho);
    }
    double tau_pspg_calc = (tau_supg / rho) * 0.02;

    // 🌟 4. 终极护甲：PSPG 绝对下限（Floor Limiter）！
    // 无论 h 有多小，tau_pspg 绝对不能低于这个保命值，否则矩阵死无全尸！
    // 如果发现 MUMPS 还是挂，就把 1e-6 往上调到 1e-5 或 1e-4。
    double min_tau_pspg = 1e-13;

    // 强行取大值！
    double tau_pspg = std::max(tau_pspg_calc, min_tau_pspg);
    // 尝试 0.5, 0.2, 或者 0.1。值越小，中心峰值速度越高（越接近理论抛物线）
    double tune_supg = 0.5e-1;
    tau_supg *= tune_supg;
    /// 强行压制PSPG

//    tbox::pout<<tau_pspg<<" aaa"<<e_Temperature<<endl;
//    TBOX_ASSERT(tau_pspg>0);
    for(int i = 0; i < n_nodes; ++i){
      double Ni = bas_val[l][i];
      double dNi_dx = bas_grad[l][i][0];
      double dNi_dy = bas_grad[l][i][1];
      double dNi_dz = bas_grad[l][i][2];
      // 测试函数沿流线方向的偏导 (用于 SUPG)
      double U_dot_gradNi = u_k * dNi_dx + v_k * dNi_dy + w_k * dNi_dz;
      for (int j = 0; j < n_nodes; ++j){
        double Nj = bas_val[l][j];
        double dNj_dx = bas_grad[l][j][0];
        double dNj_dy = bas_grad[l][j][1];
        double dNj_dz = bas_grad[l][j][2];
        // 试探函数沿流线方向的偏导
        double U_dot_gradNj = u_k * dNj_dx + v_k * dNj_dy + w_k * dNj_dz;
        /// 单元自由度映射
        int row_u = i * 4 + 0; int col_u = j * 4 + 0;
        int row_v = i * 4 + 1; int col_v = j * 4 + 1;
        int row_w = i * 4 + 2; int col_w = j * 4 + 2;
        int row_p = i * 4 + 3; int col_p = j * 4 + 3;
        /// 粘性项
        double diff = JxW * mu_value * (dNi_dx * dNj_dx + dNi_dy * dNj_dy + dNi_dz * dNj_dz);

        double G_x = -JxW * dNi_dx * Nj;
        double G_y = -JxW * dNi_dy * Nj;
        double G_z = -JxW * dNi_dz * Nj;

        // 散度乘以 -1 (与 G 形成对称结构)
        double D_x = -JxW * Ni * dNj_dx;
        double D_y = -JxW * Ni * dNj_dy;
        double D_z = -JxW * Ni * dNj_dz;

        /// 标准 Galerkin 牛顿非线性对流雅可比
        /// 数学形式: \rho * Ni * (U^k \cdot \nabla \delta u) + \rho * Ni * (\delta u \cdot \nabla U^k)

        double rho_Ni = JxW * rho * Ni;

        // 主对角块 (U^k 传输 \delta u，外加 \delta u 自身方向的梯度扰动)
        double conv_uu = rho_Ni * (U_dot_gradNj + Nj * du_dx);
        double conv_vv = rho_Ni * (U_dot_gradNj + Nj * dv_dy);
        double conv_ww = rho_Ni * (U_dot_gradNj + Nj * dw_dz);
#if PICARD
        conv_uu = rho_Ni * U_dot_gradNj;
        conv_vv = rho_Ni * U_dot_gradNj;
        conv_ww = rho_Ni * U_dot_gradNj;
#endif

        // 非对角交叉块 (\delta u 在其他方向上的梯度扰动)
        double conv_uv = rho_Ni * Nj * du_dy;  // v_j 扰动 u 方程
        double conv_uw = rho_Ni * Nj * du_dz;  // w_j 扰动 u 方程

        double conv_vu = rho_Ni * Nj * dv_dx;  // u_j 扰动 v 方程
        double conv_vw = rho_Ni * Nj * dv_dz;  // w_j 扰动 v 方程

        double conv_wu = rho_Ni * Nj * dw_dx;  // u_j 扰动 w 方程
        double conv_wv = rho_Ni * Nj * dw_dy;  // v_j 扰动 w 方程
#if PICARD
        conv_uv = 0.0;
        conv_uw = 0.0;
        conv_vu = 0.0;
        conv_vw = 0.0;
        conv_wu = 0.0;
        conv_wv = 0.0;
#endif
        /// ----------------------------------------------------
        /// SUPG & PSPG 稳定化项
        /// ----------------------------------------------------
        // PSPG 压力-压力项
        double pspg = -JxW * tau_pspg * (dNi_dx * dNj_dx + dNi_dy * dNj_dy + dNi_dz * dNj_dz);

        // SUPG 对流稳定项 (沿流线方向增加人工扩散)
        // \int \tau_{supg} (U^k \cdot \nabla N_i) * \rho * (U^k \cdot \nabla N_j)
        double supg_conv_base = JxW * tau_supg * rho * U_dot_gradNi;

        // SUPG 作用于牛顿雅可比的两部分偏导
        /// SUPG项只有对角位置有耦合，才有U_dot_gradNj
        double supg_uu = supg_conv_base * (U_dot_gradNj + Nj * du_dx);
        double supg_vv = supg_conv_base * (U_dot_gradNj + Nj * dv_dy);
        double supg_ww = supg_conv_base * (U_dot_gradNj + Nj * dw_dz);

        double supg_uv = supg_conv_base * (Nj * du_dy);
        double supg_uw = supg_conv_base * (Nj * du_dz);
        double supg_vu = supg_conv_base * (Nj * dv_dx);
        double supg_vw = supg_conv_base * (Nj * dv_dz);
        double supg_wu = supg_conv_base * (Nj * dw_dx);
        double supg_wv = supg_conv_base * (Nj * dw_dy);
        #if PICARD
        supg_uu = supg_conv_base * (U_dot_gradNj );
        supg_vv = supg_conv_base * (U_dot_gradNj );
        supg_ww = supg_conv_base * (U_dot_gradNj );

        supg_uv = 0;
        supg_uw = 0;
        supg_vu = 0;
        supg_vw = 0;
        supg_wu = 0;
        supg_wv = 0;
        #endif

        // SUPG 压力梯度项
        // \int \tau_{supg} (U^k \cdot \nabla N_i) * \nabla N_j
        double supg_p_x = JxW * tau_supg * U_dot_gradNi * dNj_dx;
        double supg_p_y = JxW * tau_supg * U_dot_gradNi * dNj_dy;
        double supg_p_z = JxW * tau_supg * U_dot_gradNi * dNj_dz;

        // 补全 PSPG 对 U, V, W 的交叉偏导数
        /// gemini让加的不知道对不对
        double pspg_u = -JxW * tau_pspg * rho * (dNi_dx * (Nj * du_dx + U_dot_gradNj) + dNi_dy * (Nj * dv_dx) + dNi_dz * (Nj * dw_dx));
        double pspg_v = -JxW * tau_pspg * rho * (dNi_dx * (Nj * du_dy) + dNi_dy * (Nj * dv_dy + U_dot_gradNj) + dNi_dz * (Nj * dw_dy));
        double pspg_w = -JxW * tau_pspg * rho * (dNi_dx * (Nj * du_dz) + dNi_dy * (Nj * dv_dz) + dNi_dz * (Nj * dw_dz + U_dot_gradNj));
#if PICARD
pspg_u = 0;
pspg_v = 0;
pspg_w = 0;
#endif
        /// ----------------------------------------------------
        /// 将所有项累加到 16x16 局部矩阵
        /// ----------------------------------------------------

        // --- 动量方程 X ---
        (*ele_mat)(row_u, col_u) += diff + conv_uu + supg_uu;
        (*ele_mat)(row_u, col_v) += conv_uv + supg_uv;
        (*ele_mat)(row_u, col_w) += conv_uw + supg_uw;
        (*ele_mat)(row_u, col_p) += G_x + supg_p_x;

        // --- 动量方程 Y ---
        (*ele_mat)(row_v, col_u) += conv_vu + supg_vu;
        (*ele_mat)(row_v, col_v) += diff + conv_vv + supg_vv;
        (*ele_mat)(row_v, col_w) += conv_vw + supg_vw;
        (*ele_mat)(row_v, col_p) += G_y + supg_p_y;

        // --- 动量方程 Z ---
        (*ele_mat)(row_w, col_u) += conv_wu + supg_wu;
        (*ele_mat)(row_w, col_v) += conv_wv + supg_wv;
        (*ele_mat)(row_w, col_w) += diff + conv_ww + supg_ww;
        (*ele_mat)(row_w, col_p) += G_z + supg_p_z;

        // --- 连续性方程 ---
//        (*ele_mat)(row_p, col_u) += D_x;
//        (*ele_mat)(row_p, col_v) += D_y;
//        (*ele_mat)(row_p, col_w) += D_z;
//        (*ele_mat)(row_p, col_p) += pspg;
        // --- 连续性方程 --- (替换你原来的这部分)
        (*ele_mat)(row_p, col_u) += D_x + pspg_u; // 加上了 pspg_u
        (*ele_mat)(row_p, col_v) += D_y + pspg_v; // 加上了 pspg_v
        (*ele_mat)(row_p, col_w) += D_z + pspg_w; // 加上了 pspg_w
        (*ele_mat)(row_p, col_p) += pspg;

      }

    }
  }
}

void LinearPrism::buildFluidResidualElementVector(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    int entity_id, tbox::Array<hier::DoubleVector<NDIM> > U_val,
    tbox::Array<double> P_val, tbox::Array<double> T_val) {
  /// Done on 2026-04-15
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager
      = IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator
      = integrator_manager->getIntegrator("LinearPrism");
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager
      = ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func
      = shape_manager->getShapeFunction("LinearPrism");
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();
  tbox::Pointer<Material> material
      = material_manager->getMaterial(GET_USER_MAT(entity_id));
  double e_Temperature = 0;
  for(int i = 0; i < 6; i++) {
    e_Temperature += T_val[i] / 6.0;
  }
  double mu_value = 1.002e-3;
  double rho = material->getDensity(e_Temperature);
  int n_nodes = shape_func->getNumberOfDof(); // 6
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  int n_dof_total = n_nodes * (NDIM + 1);     // 24

  tbox::Array<double> weight = integrator->getQuadratureWeights();

  double volume = integrator->getElementVolume();
  /// 与四面体不同，这个地方只需要它的局部积分点信息
  tbox::Array<hier::DoubleVector<NDIM> > local_quad_pnt =
      integrator->getLocalQuadraturePoints();
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, local_quad_pnt);
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, local_quad_pnt);
  for (int l = 0; l < num_quad_pnts; ++l){
    /// 动态计算雅可比行列式 |J|
    double jac = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
    /// 真实的物理体积微元 JxW
    double JxW = jac * weight[l];
    /// 3. Shakib-Hughes 动态各向异性特征长度与 PSPG/SUPG 稳定化
    double u_k = 0.0, v_k = 0.0, w_k = 0.0, p_k = 0.0;
    double du_dx = 0.0, du_dy = 0.0, du_dz = 0.0;
    double dv_dx = 0.0, dv_dy = 0.0, dv_dz = 0.0;
    double dw_dx = 0.0, dw_dy = 0.0, dw_dz = 0.0;
    double dp_dx = 0.0, dp_dy = 0.0, dp_dz = 0.0;

    for (int m = 0; m < n_nodes; ++m) {
      double N_m = bas_val[l][m];
      double dNm_dx = bas_grad[l][m][0];
      double dNm_dy = bas_grad[l][m][1];
      double dNm_dz = bas_grad[l][m][2];

      u_k += N_m * U_val[m][0];
      v_k += N_m * U_val[m][1];
      w_k += N_m * U_val[m][2];
      p_k += N_m * P_val[m];

      du_dx += dNm_dx * U_val[m][0]; du_dy += dNm_dy * U_val[m][0]; du_dz += dNm_dz * U_val[m][0];
      dv_dx += dNm_dx * U_val[m][1]; dv_dy += dNm_dy * U_val[m][1]; dv_dz += dNm_dz * U_val[m][1];
      dw_dx += dNm_dx * U_val[m][2]; dw_dy += dNm_dy * U_val[m][2]; dw_dz += dNm_dz * U_val[m][2];
      dp_dx += dNm_dx * P_val[m]; dp_dy += dNm_dy * P_val[m]; dp_dz += dNm_dz * P_val[m];

    }

    double U_norm = sqrt(u_k * u_k + v_k * v_k + w_k * w_k);
    /// ====================================================================
    /// 第 2 步：拿着完整的流速 (u_k, v_k, w_k)，开一个新循环计算特征尺度参数
    /// ====================================================================
    double u_grad_N_sum = 0.0;
    double grad_N_squared_sum = 0.0;

    for (int m = 0; m < n_nodes; ++m) {
      double dNm_dx = bas_grad[l][m][0];
      double dNm_dy = bas_grad[l][m][1];
      double dNm_dz = bas_grad[l][m][2];

      // 核心：完整的流速向量 与 当前节点形函数梯度的内积 的绝对值
      u_grad_N_sum += fabs(u_k * dNm_dx + v_k * dNm_dy + w_k * dNm_dz);

      // 顺便把 PSPG 需要的纯几何对角线也累加了
      grad_N_squared_sum += (dNm_dx * dNm_dx + dNm_dy * dNm_dy + dNm_dz * dNm_dz);
    }



    double h_supg = 0.0;
    double h_pspg = 0.0;
    /// 对于压力连续性方程，它不具备流动方向性偏好，
    /// 使用基于形函数几何梯度的调和平均，自动捕捉三棱柱的“最短边（厚度）”作为防过度耗散的基准尺度
    if (grad_N_squared_sum > 1e-15) {
      h_pspg = 2.0 / sqrt(grad_N_squared_sum);
    } else {
      h_pspg = pow(volume * jac, 1.0 / 3.0); // 极小概率退化保护
    }

    /// PSPG 稳定化参数
//    double tau_pspg = (h_pspg * h_pspg) / (12.0 * mu_value);

    /// (未来用于 N-S Jacobian，这里计算 tau_supg)
    double tau_supg = 0.0;
    if (U_norm > 1e-12 && u_grad_N_sum > 1e-12) {
      h_supg = 2.0 * U_norm / u_grad_N_sum;
      tau_supg = 1.0 / sqrt(pow(2.0 * U_norm / h_supg, 2.0) + pow(4.0 * mu_value / (rho * h_supg * h_supg), 2.0));
    }
    else {
      // 纯扩散极限下，复用 h_pspg 作为特征尺度，避免截断为 0
      tau_supg = (h_pspg * h_pspg) / (4.0 * mu_value / rho);
    }
    // 尝试 0.5, 0.2, 或者 0.1。值越小，中心峰值速度越高（越接近理论抛物线）
    double tau_pspg_calc = (tau_supg / rho) * 0.02;

    // 4. 终极护甲：PSPG 绝对下限（Floor Limiter）！
    // 无论 h 有多小，tau_pspg 绝对不能低于这个保命值，否则矩阵死无全尸！
    // 如果发现 MUMPS 还是挂，就把 1e-6 往上调到 1e-5 或 1e-4。
    double min_tau_pspg = 1e-13;

    // 强行取大值！
    double tau_pspg = std::max(tau_pspg_calc, min_tau_pspg);
    double tune_supg = 0.5e-1;
    tau_supg *= tune_supg;
    /// 强行压制PSPG

    /// 结束SUPG和PSPG参数的确定
    /// =========================================================
    /// 计算强形式动量残差
    /// 注意: 对于线性四面体 (P1)，速度的二阶导数拉普拉斯项为 0，所以粘性残差消失。
    /// 注意: 在工程实践中，为了避免极其昂贵的二阶导数计算，
    /// 对于所有一次单元（Tet4, Prism6, Hex8），我们都在 SUPG/PSPG 的
    /// 强残差中忽略了粘性拉普拉斯项 (-\mu \Delta u)。
    /// 故 Mom_R \approx \rho (u \cdot \nabla u) + \nabla p
    /// Mom_R = \rho (u \cdot \nabla u) + \nabla p
    /// =========================================================
    double Mom_Rx = rho * (u_k * du_dx + v_k * du_dy + w_k * du_dz) + dp_dx;
    double Mom_Ry = rho * (u_k * dv_dx + v_k * dv_dy + w_k * dv_dz) + dp_dy;
    double Mom_Rz = rho * (u_k * dw_dx + v_k * dw_dy + w_k * dw_dz) + dp_dz;

//    /// 修改为：(仅用于 SUPG 和 PSPG 的残差计算，主 Galerkin 残差绝对不要动！)
//    /// 防止SUPG也出现psedo扩散
//    Mom_Rx = rho * (u_k * du_dx + v_k * du_dy + w_k * du_dz);
//    Mom_Ry = rho * (u_k * dv_dx + v_k * dv_dy + w_k * dv_dz);
//    Mom_Rz = rho * (u_k * dw_dx + v_k * dw_dy + w_k * dw_dz);

    /// 连续性强残差: div(u)
    double div_u = du_dx + dv_dy + dw_dz;
    /// 组装积分
    for (int i = 0; i < n_nodes; ++i){
      double Ni = bas_val[l][i];
      double dNi_dx = bas_grad[l][i][0];
      double dNi_dy = bas_grad[l][i][1];
      double dNi_dz = bas_grad[l][i][2];

      int row_u = i * 4 + 0;
      int row_v = i * 4 + 1;
      int row_w = i * 4 + 2;
      int row_p = i * 4 + 3;

      /// ----------------------------------------------------
      /// 1. 动量方程 (X, Y, Z) - 取【负残差】 (-R_u)
      /// ----------------------------------------------------
      // a. 纯粘性扩散积分: -\int \mu \nabla u : \nabla v
      double diff_u = -JxW * mu_value * (dNi_dx * du_dx + dNi_dy * du_dy + dNi_dz * du_dz);
      double diff_v = -JxW * mu_value * (dNi_dx * dv_dx + dNi_dy * dv_dy + dNi_dz * dv_dz);
      double diff_w = -JxW * mu_value * (dNi_dx * dw_dx + dNi_dy * dw_dy + dNi_dz * dw_dz);

      // b. 纯对流平流积分: -\int \rho (u \cdot \nabla u) v
      double conv_u = -JxW * Ni * rho * (u_k * du_dx + v_k * du_dy + w_k * du_dz);
      double conv_v = -JxW * Ni * rho * (u_k * dv_dx + v_k * dv_dy + w_k * dv_dz);
      double conv_w = -JxW * Ni * rho * (u_k * dw_dx + v_k * dw_dy + w_k * dw_dz);

      // c. 压力梯度积分: +\int p (\nabla \cdot v)  (弱形式分部积分翻转符号)
      double gradp_u = JxW * p_k * dNi_dx;
      double gradp_v = JxW * p_k * dNi_dy;
      double gradp_w = JxW * p_k * dNi_dz;

      // d. SUPG 稳定化负残差: -\int \tau_{supg} (u \cdot \nabla N_i) * (Mom_R)
      double U_dot_gradNi = u_k * dNi_dx + v_k * dNi_dy + w_k * dNi_dz;
      double supg_u = -JxW * tau_supg * U_dot_gradNi * Mom_Rx;
      double supg_v = -JxW * tau_supg * U_dot_gradNi * Mom_Ry;
      double supg_w = -JxW * tau_supg * U_dot_gradNi * Mom_Rz;

      /// 累加动量右端项
      (*ele_vec)[row_u] += diff_u + conv_u + gradp_u + supg_u;
      (*ele_vec)[row_v] += diff_v + conv_v + gradp_v + supg_v;
      (*ele_vec)[row_w] += diff_w + conv_w + gradp_w + supg_w;

      /// ----------------------------------------------------
      /// 2. 连续性方程 (压力 P) - 取【正残差】 (+R_p)
      /// ----------------------------------------------------
      // a. 散度积分: +\int q (\nabla \cdot u)
      double div_term = JxW * Ni * div_u;

      // b. PSPG 稳定化正残差: +\int \tau_{pspg} (\nabla q \cdot Mom_R)
      double pspg_term = JxW * tau_pspg * (dNi_dx * Mom_Rx + dNi_dy * Mom_Ry + dNi_dz * Mom_Rz);

      /// 累加连续性方程右端项
      (*ele_vec)[row_p] += div_term + pspg_term;

    }
  }


}

double LinearPrism::calcDynamicDetJ(
    const tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& local_pnt) {

  double xi   = local_pnt[0];
  double eta  = local_pnt[1];
  double zeta = local_pnt[2];

  // 1. 获取局部偏导数
  double dN_local[6][3];
  dN_local[0][0] = -(1.0 - zeta) / 2.0; dN_local[0][1] = -(1.0 - zeta) / 2.0; dN_local[0][2] = -(1.0 - xi - eta) / 2.0;
  dN_local[1][0] =  (1.0 - zeta) / 2.0; dN_local[1][1] =  0.0;                dN_local[1][2] = -xi / 2.0;
  dN_local[2][0] =  0.0;                dN_local[2][1] =  (1.0 - zeta) / 2.0; dN_local[2][2] = -eta / 2.0;
  dN_local[3][0] = -(1.0 + zeta) / 2.0; dN_local[3][1] = -(1.0 + zeta) / 2.0; dN_local[3][2] =  (1.0 - xi - eta) / 2.0;
  dN_local[4][0] =  (1.0 + zeta) / 2.0; dN_local[4][1] =  0.0;                dN_local[4][2] =  xi / 2.0;
  dN_local[5][0] =  0.0;                dN_local[5][1] =  (1.0 + zeta) / 2.0; dN_local[5][2] =  eta / 2.0;

  // 2. 组装转置雅可比矩阵
  double J[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 6; ++k) {
        J[i][j] += dN_local[k][i] * real_vertex[k][j];
      }
    }
  }

  // 3. 计算并返回绝对值 |detJ|
  double detJ = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
              - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
              + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

  return std::abs(detJ);
}

const int LinearPrism::getProblemDim() { return NDIM; }

//获取实体自由度
tbox::Array<int> LinearPrism::getNumberOfDofOnEntity() {
  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearPrism");
  return shape_func->getNumberOfDofOnEntity();
}

