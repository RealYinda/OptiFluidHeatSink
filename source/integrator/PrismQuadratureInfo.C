/// 文件名:     PrismQuadratureInfo.C


#include "PrismQuadratureInfo.h"


/*****************************************************************************
 * 构造函数.
 *****************************************************************************/
PrismQuadratureInfo::PrismQuadratureInfo(int algebric_accuracy) {
    // 三棱柱参考单元体积: 底面积(0.5) * 高度[-1, 1](2.0) = 1.0
    d_volume = 1.0;

    if (algebric_accuracy == 1) {
        /// 若积分精度为 1 则设置积分信息为 1 点积分 (中心点):
        d_number_quad_points = 1;
        d_algebric_accuracy = algebric_accuracy;
        d_quad_points.resizeArray(1);
        d_quad_weights.resizeArray(1);

        for (int i = 0; i < d_number_quad_points; ++i) {
            d_quad_weights[i] = 1.0; // 权重和等于参考体积 1.0
        }

        // 三棱柱中心点坐标
        d_quad_points[0][0] = 1.0 / 3.0;
        d_quad_points[0][1] = 1.0 / 3.0;
        d_quad_points[0][2] = 0.0;
    }

    if (algebric_accuracy >= 2) {
        /// 若积分精度为 2 (及以上) 则设置积分信息为 6 点积分:
        d_number_quad_points = 6;
        d_algebric_accuracy = algebric_accuracy;
        d_quad_points.resizeArray(6);
        d_quad_weights.resizeArray(6);

        for (int i = 0; i < d_number_quad_points; ++i) {
            d_quad_weights[i] = 1.0 / 6.0; // 权重和为 6 * (1/6) = 1.0
        }

        double sqrt3_inv = 0.57735026918962576451; // 1 / sqrt(3)

        // 下半层积分点 (zeta = -1 / sqrt(3))
        d_quad_points[0][0] = 1.0 / 6.0;
        d_quad_points[0][1] = 1.0 / 6.0;
        d_quad_points[0][2] = -sqrt3_inv;

        d_quad_points[1][0] = 2.0 / 3.0;
        d_quad_points[1][1] = 1.0 / 6.0;
        d_quad_points[1][2] = -sqrt3_inv;

        d_quad_points[2][0] = 1.0 / 6.0;
        d_quad_points[2][1] = 2.0 / 3.0;
        d_quad_points[2][2] = -sqrt3_inv;

        // 上半层积分点 (zeta = +1 / sqrt(3))
        d_quad_points[3][0] = 1.0 / 6.0;
        d_quad_points[3][1] = 1.0 / 6.0;
        d_quad_points[3][2] = sqrt3_inv;

        d_quad_points[4][0] = 2.0 / 3.0;
        d_quad_points[4][1] = 1.0 / 6.0;
        d_quad_points[4][2] = sqrt3_inv;

        d_quad_points[5][0] = 1.0 / 6.0;
        d_quad_points[5][1] = 2.0 / 3.0;
        d_quad_points[5][2] = sqrt3_inv;
    }
}

/*****************************************************************************
 * 析构函数.
 *****************************************************************************/
PrismQuadratureInfo::~PrismQuadratureInfo() {}

/*****************************************************************************
 * 获取积分点.
 *****************************************************************************/
const tbox::Array<hier::DoubleVector<NDIM> >&
PrismQuadratureInfo::getQuadraturePoints() {
    return d_quad_points;
}

/*****************************************************************************
 * 获取积分点权重.
 *****************************************************************************/
const tbox::Array<double>& PrismQuadratureInfo::getQuadratureWeights() {
    return d_quad_weights;
}

/*****************************************************************************
 * 获取积分点数目.
 *****************************************************************************/
int PrismQuadratureInfo::getNumberOfQuadraturePoints() {
    return d_number_quad_points;
}

/*****************************************************************************
 * 获取积分单元体积.
 *****************************************************************************/
double PrismQuadratureInfo::getElementVolume() { return d_volume; }

