/// 文件名:     PrismIntegrator.h

#include "DoubleVector.h"
#include "PrismShapeFunction.h"
#include "PrismIntegrator.h"

/****************************************************************************
 * 构造函数.
 ****************************************************************************/
PrismIntegrator::PrismIntegrator(int algebric_accuracy,
                                 const string &name)
    : BaseIntegrator<NDIM>(name) {
    d_num_vertex = 6;
    d_template_vertex.resizeArray(d_num_vertex);

    // 初始化三棱柱模板单元结点 (严格匹配之前讨论过的局部参考系)
    // 底面 (zeta = -1)
    d_template_vertex[0][0] = 0.0;
    d_template_vertex[0][1] = 0.0;
    d_template_vertex[0][2] = -1.0;
    d_template_vertex[1][0] = 1.0;
    d_template_vertex[1][1] = 0.0;
    d_template_vertex[1][2] = -1.0;
    d_template_vertex[2][0] = 0.0;
    d_template_vertex[2][1] = 1.0;
    d_template_vertex[2][2] = -1.0;
    // 顶面 (zeta = 1)
    d_template_vertex[3][0] = 0.0;
    d_template_vertex[3][1] = 0.0;
    d_template_vertex[3][2] = 1.0;
    d_template_vertex[4][0] = 1.0;
    d_template_vertex[4][1] = 0.0;
    d_template_vertex[4][2] = 1.0;
    d_template_vertex[5][0] = 0.0;
    d_template_vertex[5][1] = 1.0;
    d_template_vertex[5][2] = 1.0;

    /// 坐标变换.
    d_coord_transform = new PrismCoordTran(d_template_vertex);
    /// 积分信息.
    d_quad_info = new PrismQuadratureInfo(algebric_accuracy);
}
/****************************************************************************
 * 析构函数.
 ****************************************************************************/
PrismIntegrator::~PrismIntegrator() {}

/****************************************************************************
 * 获取积分点. (注：返回的是转换到物理空间的全局积分点)
* 从纯数学矩阵组装的角度看，ShapeFunction 本质上只需要局部坐标(xi, eta, zeta)
 * 即可计算出形函数的值，将局部坐标转为全局坐标再传给 ShapeFunction 似乎是冗余的
 * (甚至会导致非线性单元如 Prism 在 ShapeFunction 内部被迫做牛顿迭代反解局部坐标)。
 *
 * 但是，在 JAUMIN 这样的工业级大型通用有限元框架中：
 * 1. Integrator 接口的定义必须具有高度的通用性，它不仅仅服务于矩阵组装(Assembly)，
 * 还要服务于流固耦合插值、粒子追踪、后处理等各种物理空间(Physical Space)操作。
 * 2. 框架的顶层契约规定了 ShapeFunction::value() 必须接收真实的全局物理坐标 (x,y,z)。
 *
 * 因此，这里必须调用 local2Global() 将 QuadratureInfo 提供的局部高斯点
 * 显式地转换为全局坐标(g_pnt)再抛出。这是一种典型的“以牺牲少量底层性能，换取整个
 * 框架 API 绝对统一和向后兼容性”的架构妥协。
 ****************************************************************************/
/// 这里需要注意的是，我传给value和gradient的参数都是局部积分点信息
/// !!!!!!!
/// 如果没有特殊需求，在三棱柱单元里这个函数不太会用到
tbox::Array<hier::DoubleVector<NDIM> >
PrismIntegrator::getQuadraturePoints(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {
  tbox::Array<hier::DoubleVector<NDIM> > l_pnt =
      d_quad_info->getQuadraturePoints();
  int num_pnt = l_pnt.getSize();
  tbox::Array<hier::DoubleVector<NDIM> > g_pnt(num_pnt);
  for (int i = 0; i < num_pnt; ++i) {
    d_coord_transform->local2Global(real_vertex, l_pnt[i], g_pnt[i]);
  }
  return g_pnt;
}

/****************************************************************************
 * 获取积分点权重.
 ****************************************************************************/
const tbox::Array<double> &PrismIntegrator::getQuadratureWeights() {
  return d_quad_info->getQuadratureWeights();
}

/****************************************************************************
 * 获取积分点数目.
 ****************************************************************************/
int PrismIntegrator::getNumberOfQuadraturePoints() {
  return d_quad_info->getNumberOfQuadraturePoints();
}

/****************************************************************************
 * 获取单元的面积(体积).
 ****************************************************************************/
double PrismIntegrator::getElementVolume() {
  return d_quad_info->getElementVolume();
}

/****************************************************************************
 * 获取仿射变换的矩阵行列式 (Local -> Global).
 ****************************************************************************/
double PrismIntegrator::getLocal2GlobalJacobian(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {
  return d_coord_transform->Local2GlobalJacobian(real_vertex);
}

/****************************************************************************
 * 获取仿射变换的矩阵行列式 (Global -> Local).
 ****************************************************************************/
double PrismIntegrator::getGlobal2LocalJacobian(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {
  return d_coord_transform->Global2LocalJacobian(real_vertex);
}

/****************************************************************************
 * 获取结点数目.
 ****************************************************************************/
int PrismIntegrator::getNumberOfVertex() { return d_num_vertex; }


/// 用在prism里
tbox::Array<hier::DoubleVector<NDIM> >PrismIntegrator::getLocalQuadraturePoints(){
  return d_quad_info->getQuadraturePoints();
}
