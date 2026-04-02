/// 文件名:     LinearPrism.C


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
  // TODO: 由用户实现右端项向量计算
}

void LinearPrism::buildTh_ElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val) {
  // TODO: 由用户实现瞬态热传导矩阵
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
        integrator_manager->getIntegrator("PrismIntegrator");

    /// 取出三棱柱形函数对象
    tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
        ShapeFunctionManager<NDIM>::getManager();
    tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
        shape_manager->getShapeFunction("PrismShapeFunction");

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
      double jac = integrator->getLocal2GlobalJacobian(real_vertex);
      tbox::Array<double> weight = integrator->getQuadratureWeights();

      /// 取出形函数在各个积分点的梯度值
      tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
          shape_func->gradient(real_vertex, local_quad_pnt);

      /// 根据平均温度获取当前材料的热导率 K
      double K = material->getK(e_Temperature);

}

void LinearPrism::buildTh_ElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    int entity_id, tbox::Array<double> T_val, double e_ThermalSource) {
  // TODO: 由用户实现瞬态热传导右端项
}

void LinearPrism::buildStaticTh_ElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    int entity_id, tbox::Array<double> T_val, double e_ThermalSource) {
  // TODO: 由用户实现稳态热传导右端项
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
