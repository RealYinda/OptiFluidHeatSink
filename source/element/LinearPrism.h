/// 文件名:     LinearPrism.h
#ifndef include_LinearPrism
#define include_LinearPrism
#include "BaseElement.h"
#include "BaseIntegrator.h"
#include "BaseMaterial.h"
#include "BaseShapeFunction.h"
#include "BaseMaterial.h"
#include "DoubleVector.h"
#include "Vector.h"
#include "Matrix.h"

#include "IntegratorManager.h"
#include "ShapeFunctionManager.h"
#include "ElementManager.h"
#include "MaterialManager.h"
#include "Material.h"

//update
#include "Pointer.h"
#include "Array.h"
//自定义宏管理文件
#include "MacrosManager.h"
#include "JAUMIN_Macros.h"

using namespace JAUMIN;
class LinearPrism : public BaseElement<NDIM> {
public:
  /**
   * @brief 构造函数.
   *
   * @param name  单元名字.
   */
  LinearPrism(const string& name);

  /**
   * @brief 析构函数.
   */
  ~LinearPrism();
  virtual void buildStiffElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, double T_val);

  virtual void buildRecoveryMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id);

  virtual void buildMassElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, double T_val);

  virtual void buildDumpElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, double T_val);

  virtual void buildElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val);

  virtual void buildElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
      NewmarkData *d_newmark, int entity_id, tbox::Array<double> T_val,
      tbox::Array<double> Tolder_val);

  virtual void buildTh_ElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val);

  virtual void buildStaticTh_ElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val);

  virtual void buildTh_ElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
      int entity_id, tbox::Array<double> T_val, double e_ThermalSource);

  virtual void buildStaticTh_ElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
      int entity_id, tbox::Array<double> T_val, double e_ThermalSource);

  virtual void buildE_ElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val);

  virtual void buildE_ElementRHS(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
      int entity_id, tbox::Array<double> T_val);

  virtual void buildInitFluidElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat, int entity_id, tbox::Array<double> T_val, double *mu);

  virtual void buildFluidJacobianElementMatrix(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,
      int entity_id, tbox::Array<double> T_val,
      tbox::Array<hier::DoubleVector<NDIM> > U_val);

  virtual void buildFluidResidualElementVector(
      tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
      const double time, tbox::Array<double>& ele_vec,
      int entity_id, tbox::Array<hier::DoubleVector<NDIM> > U_val,
      tbox::Array<double> P_val, tbox::Array<double> T_val);

  /**
     * @brief 获取单元求解问题维数.
     */
  const int getProblemDim();

  /**
     * @brief 获取单元上每一个网格实体上自由度的数目.
     */
  tbox::Array<int> getNumberOfDofOnEntity();
/// 私有类，用于三棱柱单元的计算
private:
  /**
     * @brief 计算三棱柱在某一个局部等参点上的动态雅可比行列式 |detJ|
     */
  double calcDynamicDetJ(
      const tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const hier::DoubleVector<NDIM>& local_pnt);

};

#endif
