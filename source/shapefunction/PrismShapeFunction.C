/// 文件名:     PrismShapeFunction.h
/// Yin-Da Wang
/// 描述　:     三棱柱有限元形函数类
#include "PrismShapeFunction.h"
/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
PrismShapeFunction::PrismShapeFunction(const string& name)
  : BaseShapeFunction<NDIM>(name) {
  d_num_dof = 6;
  d_num_dof_on_entity.resizeArray(4);
  // 严格修正：每个顶点 (Vertex) 上仅附着 1 个自由度
  d_num_dof_on_entity[0] = 1;
  d_num_dof_on_entity[1] = 0;
  d_num_dof_on_entity[2] = 0;
  d_num_dof_on_entity[3] = 0;
}
/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
PrismShapeFunction::~PrismShapeFunction() {}

/*****************************************************************************
 * 获取某一个基函数在某一点的值.
 *****************************************************************************/
double PrismShapeFunction::value(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<double> all_vals = value(real_vertex, pnt);
  return all_vals[j];
}
/*****************************************************************************
 * 获取某一个基函数在某些点的值.
 *****************************************************************************/
tbox::Array<double> PrismShapeFunction::value(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<double> tmp_val(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tmp_val[i] = value(j, real_vertex, pnts[i]);
  }
  return tmp_val;
}

/*****************************************************************************
 * 获取基函数在某一点的值.
 *****************************************************************************/
tbox::Array<double> PrismShapeFunction::value(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<double> tmp_val(6);

  // 直接提取局部坐标
  double xi   = pnt[0];
  double eta  = pnt[1];
  double zeta = pnt[2];

  // 代入张量积形函数公式 (极其高效的纯代数计算)
  // 底面节点 (zeta = -1)
  tmp_val[0] = xi * (1.0 - zeta) / 2.0;
  tmp_val[1] = eta * (1.0 - zeta) / 2.0;
  tmp_val[2] = (1.0 - xi - eta) * (1.0 - zeta) / 2.0;

  // 顶面节点 (zeta = 1)
  tmp_val[3] = xi * (1.0 + zeta) / 2.0;
  tmp_val[4] = eta * (1.0 + zeta) / 2.0;
  tmp_val[5] = (1.0 - xi - eta) * (1.0 + zeta) / 2.0;

  return tmp_val;
}

/*****************************************************************************
 * 获取基函数在某些点的值.
 *****************************************************************************/
tbox::Array<tbox::Array<double> > PrismShapeFunction::value(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<tbox::Array<double> > tmp_val(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tbox::Array<double> val_i = value(real_vertex, pnts[i]);
    tmp_val[i] = val_i;
  }
  return tmp_val;
}

/*****************************************************************************
 * 获取某一个基函数在某一点的梯度值.
 *****************************************************************************/
tbox::Array<double> PrismShapeFunction::gradient(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<tbox::Array<double> > all_grads = gradient(real_vertex, pnt);
  return all_grads[j];
}

/*****************************************************************************
 * 获取某一个基函数在某些点的梯度值.
 *****************************************************************************/
tbox::Array<tbox::Array<double> > PrismShapeFunction::gradient(
    int j, tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<tbox::Array<double> > tmp_grad(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tbox::Array<double> grad = gradient(j, real_vertex, pnts[i]);
    tmp_grad[i] = grad;
  }
  return tmp_grad;
}

/*****************************************************************************
 * 获取基函数在某一点的梯度值 (动态计算雅可比矩阵映射).
 *****************************************************************************/
tbox::Array<tbox::Array<double> > PrismShapeFunction::gradient(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const hier::DoubleVector<NDIM>& pnt) {
  tbox::Array<tbox::Array<double> > tmp_grad(6);

  double xi   = pnt[0];
  double eta  = pnt[1];
  double zeta = pnt[2];

  // 1. 计算对局部坐标 (xi, eta, zeta) 的偏导数
  double dN_local[6][3];

  dN_local[0][0] =  (1.0 - zeta) / 2.0;
  dN_local[0][1] =  0.0;
  dN_local[0][2] = -xi / 2.0;

  dN_local[1][0] =  0.0;
  dN_local[1][1] =  (1.0 - zeta) / 2.0;
  dN_local[1][2] = -eta / 2.0;

  dN_local[2][0] = -(1.0 - zeta) / 2.0;
  dN_local[2][1] = -(1.0 - zeta) / 2.0;
  dN_local[2][2] = -(1.0 - xi - eta) / 2.0;

  dN_local[3][0] =  (1.0 + zeta) / 2.0;
  dN_local[3][1] =  0.0;
  dN_local[3][2] =  xi / 2.0;

  dN_local[4][0] =  0.0;
  dN_local[4][1] =  (1.0 + zeta) / 2.0;
  dN_local[4][2] =  eta / 2.0;

  dN_local[5][0] = -(1.0 + zeta) / 2.0;
  dN_local[5][1] = -(1.0 + zeta) / 2.0;
  dN_local[5][2] =  (1.0 - xi - eta) / 2.0;

  // 2. 组装雅可比矩阵 J = dx_i / dxi_j
  double J[3][3] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 6; ++k) {
        J[i][j] += dN_local[k][i] * real_vertex[k][j];
      }
    }
  }

  // 3. 计算行列式及伴随/逆矩阵
  double detJ = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1])
      - J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0])
      + J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);

  double invJ[3][3];
  invJ[0][0] =  (J[1][1]*J[2][2] - J[1][2]*J[2][1]) / detJ;
  invJ[0][1] = -(J[0][1]*J[2][2] - J[0][2]*J[2][1]) / detJ;
  invJ[0][2] =  (J[0][1]*J[1][2] - J[0][2]*J[1][1]) / detJ;
  invJ[1][0] = -(J[1][0]*J[2][2] - J[1][2]*J[2][0]) / detJ;
  invJ[1][1] =  (J[0][0]*J[2][2] - J[0][2]*J[2][0]) / detJ;
  invJ[1][2] = -(J[0][0]*J[1][2] - J[0][2]*J[1][0]) / detJ;
  invJ[2][0] =  (J[1][0]*J[2][1] - J[1][1]*J[2][0]) / detJ;
  invJ[2][1] = -(J[0][0]*J[2][1] - J[0][1]*J[2][0]) / detJ;
  invJ[2][2] =  (J[0][0]*J[1][1] - J[0][1]*J[1][0]) / detJ;

  // 4. 将局部偏导数通过雅可比逆矩阵映射到全局物理空间
  for (int i = 0; i < 6; ++i) {
    tmp_grad[i].resizeArray(3);
    tmp_grad[i][0] = dN_local[i][0]*invJ[0][0] + dN_local[i][1]*invJ[1][0] + dN_local[i][2]*invJ[2][0];
    tmp_grad[i][1] = dN_local[i][0]*invJ[0][1] + dN_local[i][1]*invJ[1][1] + dN_local[i][2]*invJ[2][1];
    tmp_grad[i][2] = dN_local[i][0]*invJ[0][2] + dN_local[i][1]*invJ[1][2] + dN_local[i][2]*invJ[2][2];
  }

  return tmp_grad;
}

/*****************************************************************************
 * 获取基函数在某些点的梯度值.
 *****************************************************************************/
tbox::Array<tbox::Array<tbox::Array<double> > >
PrismShapeFunction::gradient(
    tbox::Array<hier::DoubleVector<NDIM> >& real_vertex,
    const tbox::Array<hier::DoubleVector<NDIM> >& pnts) {
  int pnt_size = pnts.getSize();
  tbox::Array<tbox::Array<tbox::Array<double> > > tmp_grad(pnt_size);
  for (int i = 0; i < pnt_size; ++i) {
    tbox::Array<tbox::Array<double> > grad = gradient(real_vertex, pnts[i]);
    tmp_grad[i] = grad;
  }
  return tmp_grad;
}

/****************************************************************************
 * 获取自由度数目.
 ****************************************************************************/
int PrismShapeFunction::getNumberOfDof() { return d_num_dof; }

tbox::Array<int> PrismShapeFunction::getNumberOfDofOnEntity() {
  return d_num_dof_on_entity;
}
