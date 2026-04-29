/// 文件名:     PrismShapeFunction.h
/// Yin-Da Wang
/// 描述　:     三棱柱有限元形函数类

#ifndef included_PrismShapeFunction
#define included_PrismShapeFunction
#include "Array.h"
#include "Pointer.h"
#include "BaseShapeFunction.h"
#include "DoubleVector.h"
using namespace JAUMIN;
/*!
 * @brief 有限元三棱柱线性形函数类，为有限元计算中的单元积分计算提供形函数在积
 * 分点上的函数值和梯度值。
 *
 */

class PrismShapeFunction : public BaseShapeFunction<NDIM>{
public:
  /**
   * @brief 构造函数。
   *
   * @param name 输入参数, 字符串, 形函数的名字.
   */
  PrismShapeFunction(const string& name);

  /**
   * @brief 析构函数。
   */
  ~PrismShapeFunction();
  /**
     * @brief 取出某一个积分点上的形函数值。
     *
     * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
     * @param pnt          输入参数，坐标，积分点坐标。
     *
     * @return 双精度型数组，积分点的形函数值。
     */
  tbox::Array<double> value(tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
                            const hier::DoubleVector<NDIM>& pnt);

  /**
     * @brief 取出某一个积分点上的某一个形函数值。
     *
     * @param j            输入参数, 整型, 积分点编号.
     * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
     * @param pnt          输入参数，坐标，积分点坐标。
     *
     * @return 双精度型，积分点的形函数值。
     */
  double value(int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
               const hier::DoubleVector<NDIM>& pnt);

  /**
     * @brief 取出若干个积分点上的某一个形函数值。
     *
     * @param j            输入参数, 整型, 形函数编号.
     * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
     * @param pnts         输入参数，坐标数组，积分点坐标数组。
     *
     * @return 双精度型数组，积分点数组的形函数值。
     */
  tbox::Array<double> value(int j,
                            tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
                            const tbox::Array<hier::DoubleVector<NDIM> >& pnts);
  /**
     * @brief 取出若干个积分点上的形函数值。
     *
     * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
     * @param pnts         输入参数，坐标数组，积分点坐标数组。
     *
     * @return 双精度型，积分点数组的形函数值。
     */
  tbox::Array<tbox::Array<double> > value(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<NDIM> >& pnts);

  /**
     * @brief 取出某一个积分点上的形函数梯度。
     *
     * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
     * @param pnt          输入参数，坐标，积分点坐标。
     *
     * @return 双精度型数组，积分点的形函数梯度。
     */
  tbox::Array<tbox::Array<double> > gradient(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const hier::DoubleVector<NDIM>& pnt);

  /**
     * @brief 取出某一个积分点上的某一个形函数梯度。
     *
     * @param j            输入参数, 整型,     形函数编号.
     * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
     * @param pnt          输入参数, 坐标,     积分点坐标.
     *
     * @return 双精度型数组，积分点的形函数梯度。
     */
  tbox::Array<double> gradient(
      int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const hier::DoubleVector<NDIM>& pnt);

  /**
     * @brief 取出某若干积分点上的形函数梯度。
     * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
     * @param pnts         输入参数，坐标，    积分点坐标数组。
     *
     * @return 二维双精度型数组，积分点数组的形函数梯度。
     */
  tbox::Array<tbox::Array<tbox::Array<double> > > gradient(
      tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<NDIM> >& pnts);

  /**
     * @brief 取出某若干积分点上的某一个形函数梯度。
     *
     * @param j            输入参数, 整型,     形函数编号.
     * @param real_vertex  输入参数, 坐标数组, 单元结点坐标.
     * @param pnts         输入参数，坐标，    积分点坐标数组。
     *
     * @return 二维双精度型数组，积分点数组的形函数梯度。
     */
  virtual tbox::Array<tbox::Array<double> > gradient(
      int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
      const tbox::Array<hier::DoubleVector<NDIM> >& pnts);
  /**
       * @brief 获取单元上自由度的数目.
       *
       * @return 整型, 单元自由度的数目.
       */
  virtual int getNumberOfDof();

  tbox::Array<int> getNumberOfDofOnEntity();
private:

  tbox::Array<int> d_num_dof_on_entity;
  int d_num_dof; /**< 该形函数上自由度的数目 */
  tbox::Array<hier::DoubleVector<NDIM> >
  d_real_vertex; /**< 实际网格单元的结点坐标 */

};

#endif
