/// 文件名:     PrismIntegrator.h

#ifndef included_fem_PrismIntegrator
#define included_fem_PrismIntegrator

#include "Array.h"
#include "Pointer.h"
#include "BaseIntegrator.h"
#include "DoubleVector.h"
#include "PrismCoordTran.h"
#include "PrismQuadratureInfo.h"
using namespace JAUMIN;
/*!
 * @brief 三棱柱有限元积分器类，提供单元上有限元积分计算所需要的积分信息。该类根
 * 据用户提供的网格单元的结点信息和积分精度初始化有限元积分器并获取有限元
 * 积分计算需要的信息，其中用到了PrismQuadratureInfo，PrismCoordTran两个成员类。
 *
 * @see PrismQuadratureInfo, PrismCoordTran
 */

class PrismIntegrator : public BaseIntegrator<NDIM> {
public:
  /**
   * @brief 构造函数。
   *
   * @param algebric_accuracy 输入参数，整型，积分精度。
   * @param ele_name          输入参数，字符串，单元名称。
   *
   */
  PrismIntegrator(int algebric_accuracy, const string &ele_name);

  /**
   * @brief 析构函数。
   */
  ~PrismIntegrator();

  /**
   * @brief 获取积分点。
   *
   * @param real_vertex 输入参数, 坐标数组, 实际网格单元的结点坐标。
   * @return 坐标数组，积分点数组(全局坐标)。
   */
  virtual tbox::Array<hier::DoubleVector<NDIM> > getQuadraturePoints(
      tbox::Array<hier::DoubleVector<NDIM> > &real_vertex);

  /**
   * @brief 获取积分点权重。
   *
   * @return 双精度浮点型数组，积分点权重数组。
   */
  virtual const tbox::Array<double> &getQuadratureWeights();

  /**
   * @brief 获取积分点数目。
   *
   * @return 整型，积分点数目。
   */
  virtual int getNumberOfQuadraturePoints();

  /**
   * @brief 获取单元的面积(体积)。
   *
   * @return 双精度浮点型，单元面积(体积)。
   */
  virtual double getElementVolume();

  /**
   * @brief 获取由局部到全局仿射变换的矩阵行列式。
   *
   * @param real_vertex 输入参数, 坐标数组, 单元的结点坐标.
   *
   * @return 双精度浮点型，仿射变换的矩阵行列式。
   *
   */
  virtual double getLocal2GlobalJacobian(
      tbox::Array<hier::DoubleVector<NDIM> > &real_vertex);

  /**
   * @brief 获取由全局到局部仿射变换的矩阵行列式。
   *
   * @param real_vertex 输入参数, 坐标数组, 单元的结点坐标.
   *
   * @return 双精度浮点型，仿射变换的矩阵行列式。
   */
  virtual double getGlobal2LocalJacobian(
      tbox::Array<hier::DoubleVector<NDIM> > &real_vertex);

  /**
   * @brief 获取单元结点数目。
   *
   * @return 整型，单元结点数目。
   */
  virtual int getNumberOfVertex();

private:
  tbox::Pointer<PrismCoordTran> d_coord_transform; /**< 坐标变换 */
  tbox::Pointer<PrismQuadratureInfo> d_quad_info;  /**< 积分信息 */
  tbox::Array<hier::DoubleVector<NDIM> >
  d_template_vertex; /**< 有限元模板单元的结点 */
  int d_num_vertex;      /**< 结点数目 */
};
#endif
