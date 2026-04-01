/// 文件名:     PrismCoordTran.C
#include "PrismCoordTran.h"
#include <cmath>

/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
PrismCoordTran::PrismCoordTran(
    tbox::Array<hier::DoubleVector<NDIM> > &template_vertex) {
  d_template_vertex = template_vertex;
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
PrismCoordTran::~PrismCoordTran() {}

/*****************************************************************************
 * 把一个模板单元上的局部积分点映射到实际网格上的积分点.
 * 正向映射: 利用形函数代数方程直接求解
 *****************************************************************************/
void PrismCoordTran::local2Global(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex,
    const hier::DoubleVector<NDIM> &lp, hier::DoubleVector<NDIM> &gp) {

  double xi   = lp[0];
  double eta  = lp[1];
  double zeta = lp[2];

  double N[6];
  N[0] = xi * (1.0 - zeta) / 2.0;
  N[1] = eta * (1.0 - zeta) / 2.0;
  N[2] = (1.0 - xi - eta) * (1.0 - zeta) / 2.0;
  N[3] = xi * (1.0 + zeta) / 2.0;
  N[4] = eta * (1.0 + zeta) / 2.0;
  N[5] = (1.0 - xi - eta) * (1.0 + zeta) / 2.0;

  gp[0] = 0.0; gp[1] = 0.0; gp[2] = 0.0;
  for (int i = 0; i < 6; ++i) {
    gp[0] += N[i] * real_vertex[i][0];
    gp[1] += N[i] * real_vertex[i][1];
    gp[2] += N[i] * real_vertex[i][2];
  }
}

/*****************************************************************************
 * 把一个实际单元上的积分点映射到模板单元上的局部积分点.
 * 逆向映射: 由于三棱柱非线性，使用 Newton-Raphson 迭代求解
 *****************************************************************************/
void PrismCoordTran::global2Local(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex,
    const hier::DoubleVector<NDIM> &gp, hier::DoubleVector<NDIM> &lp) {

  double x_target = gp[0];
  double y_target = gp[1];
  double z_target = gp[2];

  // 初始猜测值设为三棱柱体心 (即 1 阶积分点位置)
  double xi   = 1.0 / 3.0;
  double eta  = 1.0 / 3.0;
  double zeta = 0.0;

  int max_iter = 15;
  double tol = 1e-12;

  for (int iter = 0; iter < max_iter; ++iter) {
    double N[6];
    N[0] = xi * (1.0 - zeta) / 2.0;
    N[1] = eta * (1.0 - zeta) / 2.0;
    N[2] = (1.0 - xi - eta) * (1.0 - zeta) / 2.0;
    N[3] = xi * (1.0 + zeta) / 2.0;
    N[4] = eta * (1.0 + zeta) / 2.0;
    N[5] = (1.0 - xi - eta) * (1.0 + zeta) / 2.0;

    double x_guess = 0.0, y_guess = 0.0, z_guess = 0.0;
    for (int i = 0; i < 6; ++i) {
      x_guess += N[i] * real_vertex[i][0];
      y_guess += N[i] * real_vertex[i][1];
      z_guess += N[i] * real_vertex[i][2];
    }

    double Fx = x_guess - x_target;
    double Fy = y_guess - y_target;
    double Fz = z_guess - z_target;

    if ((Fx*Fx + Fy*Fy + Fz*Fz) < tol) {
      break;
    }

    double dN_local[6][3];
    dN_local[0][0] =  (1.0 - zeta) / 2.0; dN_local[0][1] =  0.0;                dN_local[0][2] = -xi / 2.0;
    dN_local[1][0] =  0.0;                dN_local[1][1] =  (1.0 - zeta) / 2.0; dN_local[1][2] = -eta / 2.0;
    dN_local[2][0] = -(1.0 - zeta) / 2.0; dN_local[2][1] = -(1.0 - zeta) / 2.0; dN_local[2][2] = -(1.0 - xi - eta) / 2.0;
    dN_local[3][0] =  (1.0 + zeta) / 2.0; dN_local[3][1] =  0.0;                dN_local[3][2] =  xi / 2.0;
    dN_local[4][0] =  0.0;                dN_local[4][1] =  (1.0 + zeta) / 2.0; dN_local[4][2] =  eta / 2.0;
    dN_local[5][0] = -(1.0 + zeta) / 2.0; dN_local[5][1] = -(1.0 + zeta) / 2.0; dN_local[5][2] =  (1.0 - xi - eta) / 2.0;

    double J[3][3] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < 6; ++k) {
          J[i][j] += dN_local[k][i] * real_vertex[k][j];
        }
      }
    }

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

    xi   -= (invJ[0][0]*Fx + invJ[0][1]*Fy + invJ[0][2]*Fz);
    eta  -= (invJ[1][0]*Fx + invJ[1][1]*Fy + invJ[1][2]*Fz);
    zeta -= (invJ[2][0]*Fx + invJ[2][1]*Fy + invJ[2][2]*Fz);
  }

  lp[0] = xi;
  lp[1] = eta;
  lp[2] = zeta;
}

/*****************************************************************************
 * 从模板单元到实际网格 Jacobian 变换行列式.
 * 对于三棱柱，返回全局真实体积 / 局部模板体积 (1.0).
 * 使用三棱柱标准的三面体分割法精确计算畸变三棱柱的真实体积。
 *****************************************************************************/
double PrismCoordTran::Local2GlobalJacobian(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {

  // 模板参考三棱柱体积为 1.0
  double lvolume = 1.0;

  // 将 6 节点三棱柱拆分为 3 个四面体:
  // Tet 1: 0, 1, 2, 5
  // Tet 2: 0, 1, 5, 4
  // Tet 3: 0, 4, 5, 3
  // get_volume() 函数返回的是 6 倍的四面体体积
  double tet1 = std::abs(get_volume(real_vertex[0], real_vertex[1], real_vertex[2], real_vertex[5]));
  double tet2 = std::abs(get_volume(real_vertex[0], real_vertex[1], real_vertex[5], real_vertex[4]));
  double tet3 = std::abs(get_volume(real_vertex[0], real_vertex[4], real_vertex[5], real_vertex[3]));

  double gvolume = (tet1 + tet2 + tet3) / 6.0;

  return gvolume / lvolume;
}

/*****************************************************************************
 * 从实际网格到模板单元 Jacobian 变换行列式.
 *****************************************************************************/
double PrismCoordTran::Global2LocalJacobian(
    tbox::Array<hier::DoubleVector<NDIM> > &real_vertex) {

  double lvolume = 1.0;

  double tet1 = std::abs(get_volume(real_vertex[0], real_vertex[1], real_vertex[2], real_vertex[5]));
  double tet2 = std::abs(get_volume(real_vertex[0], real_vertex[1], real_vertex[5], real_vertex[4]));
  double tet3 = std::abs(get_volume(real_vertex[0], real_vertex[4], real_vertex[5], real_vertex[3]));

  double gvolume = (tet1 + tet2 + tet3) / 6.0;

  return lvolume / gvolume;
}

double PrismCoordTran::get_volume(const hier::DoubleVector<NDIM> &v0,
                                  const hier::DoubleVector<NDIM> &v1,
                                  const hier::DoubleVector<NDIM> &v2,
                                  const hier::DoubleVector<NDIM> &v3) {
  return ((v1[0] - v0[0]) * (v2[1] - v0[1]) * (v3[2] - v0[2]) +
          (v1[1] - v0[1]) * (v2[2] - v0[2]) * (v3[0] - v0[0]) +
          (v1[2] - v0[2]) * (v2[0] - v0[0]) * (v3[1] - v0[1]) -
          (v1[0] - v0[0]) * (v2[2] - v0[2]) * (v3[1] - v0[1]) -
          (v1[1] - v0[1]) * (v2[0] - v0[0]) * (v3[2] - v0[2]) -
          (v1[2] - v0[2]) * (v2[1] - v0[1]) * (v3[0] - v0[0]));
}
