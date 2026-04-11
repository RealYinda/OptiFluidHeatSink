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

void LinearPrism::buildStaticTh_ElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val) {
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

  /// 根据平均温度获取当前材料的热导率 K
  double K = material->getK(e_Temperature);

  for (int l = 0; l < num_quad_pnts; ++l){
    double detJ = calcDynamicDetJ(real_vertex, local_quad_pnt[l]);
    // 计算当前积分点的体积微元
    double JxW = detJ * weight[l];
    // 将刚度项累加到单元矩阵中
    for (int i = 0; i < n_dof; ++i) {
      for (int j = 0; j < n_dof; ++j) {
        (*ele_mat)(i, j) += JxW * ((bas_grad[l][i][0] * bas_grad[l][j][0] +
            bas_grad[l][i][1] * bas_grad[l][j][1] +
            bas_grad[l][i][2] * bas_grad[l][j][2]) * K);
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
      double tau_pspg = (h_pspg * h_pspg) / (12.0 * mu_value);

      /// (如果未来用于 N-S Jacobian，这里计算 tau_supg)
      double tau_supg = 0.0;
      if (U_norm > 1e-12 && u_grad_N_sum > 1e-12) {
        h_supg = 2.0 * U_norm / u_grad_N_sum;
        tau_supg = 1.0 / sqrt(pow(2.0 * U_norm / h_supg, 2.0) + pow(4.0 * mu_value / (rho * h_supg * h_supg), 2.0));
      }
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
}

void LinearPrism::buildFluidResidualElementVector(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Array<double>& ele_vec,
    int entity_id, tbox::Array<hier::DoubleVector<NDIM> > U_val,
    tbox::Array<double> P_val, tbox::Array<double> T_val) {
  // TODO: 由用户实现流体残差向量
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

