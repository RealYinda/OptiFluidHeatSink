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
