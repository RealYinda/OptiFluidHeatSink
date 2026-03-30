//
// 文件名:     LinearTet.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:39:40 2014 $
// 描述　:     四面体线弹性单元.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "LinearTet.h"

/*************************************************************************
 * 构造函数.
 ************************************************************************/
LinearTet::LinearTet(const string& name) : BaseElement<NDIM>(name) {}

/*************************************************************************
 * 析构函数.
 ************************************************************************/
LinearTet::~LinearTet() {}


/// 恢复型算法的实现
void LinearTet::buildRecoveryMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id){
  double T_val = 300.;
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  //update #3:
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");
  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");
  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);
  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);
  ele_mat->resize(n_dof * NDIM, 24);
  for (int i = 0; i < n_dof * NDIM; ++i) {
    for (int j = 0; j < 24; ++j) {
      (*ele_mat)(i, j) = 0.0;
    }
  }
  for(int quad = 0; quad <num_quad_pnts; quad++){
    double JxW = volume * jac * weight[quad];
    /// 本个积分点的坐标值
    double quad_coord[NDIM] = {0.,0.,0.};
    for(int i = 0; i < n_dof; i++){
      quad_coord[0] += real_vertex[i][0]*bas_val[quad][i];
      quad_coord[1] += real_vertex[i][1]*bas_val[quad][i];
      quad_coord[2] += real_vertex[i][2]*bas_val[quad][i];
    }
    for(int i = 0; i < n_dof; i++){
      /// 自由度i的x分量 146
      (*ele_mat)(NDIM * i + 0, 0 ) += bas_grad[quad][i][0]*quad_coord[0]*JxW;
      (*ele_mat)(NDIM * i + 0, 1 ) += bas_grad[quad][i][0]*quad_coord[1]*JxW;
      (*ele_mat)(NDIM * i + 0, 2 ) += bas_grad[quad][i][0]*quad_coord[2]*JxW;
      (*ele_mat)(NDIM * i + 0, 3 ) += bas_grad[quad][i][0]*1.*JxW;
      (*ele_mat)(NDIM * i + 0, 12 ) += bas_grad[quad][i][1]*quad_coord[0]*JxW;
      (*ele_mat)(NDIM * i + 0, 13 ) += bas_grad[quad][i][1]*quad_coord[1]*JxW;
      (*ele_mat)(NDIM * i + 0, 14 ) += bas_grad[quad][i][1]*quad_coord[2]*JxW;
      (*ele_mat)(NDIM * i + 0, 15 ) += bas_grad[quad][i][1]*1.*JxW;
      (*ele_mat)(NDIM * i + 0, 20 ) += bas_grad[quad][i][2]*quad_coord[0]*JxW;
      (*ele_mat)(NDIM * i + 0, 21 ) += bas_grad[quad][i][2]*quad_coord[1]*JxW;
      (*ele_mat)(NDIM * i + 0, 22 ) += bas_grad[quad][i][2]*quad_coord[2]*JxW;
      (*ele_mat)(NDIM * i + 0, 23 ) += bas_grad[quad][i][2]*1.*JxW;
      /// 自由度i的y分量 425
      (*ele_mat)(NDIM * i + 1, 12 ) += bas_grad[quad][i][0]*quad_coord[0]*JxW;
      (*ele_mat)(NDIM * i + 1, 13 ) += bas_grad[quad][i][0]*quad_coord[1]*JxW;
      (*ele_mat)(NDIM * i + 1, 14 ) += bas_grad[quad][i][0]*quad_coord[2]*JxW;
      (*ele_mat)(NDIM * i + 1, 15 ) += bas_grad[quad][i][0]*1.*JxW;
      (*ele_mat)(NDIM * i + 1, 4 ) += bas_grad[quad][i][1]*quad_coord[0]*JxW;
      (*ele_mat)(NDIM * i + 1, 5 ) += bas_grad[quad][i][1]*quad_coord[1]*JxW;
      (*ele_mat)(NDIM * i + 1, 6 ) += bas_grad[quad][i][1]*quad_coord[2]*JxW;
      (*ele_mat)(NDIM * i + 1, 7 ) += bas_grad[quad][i][1]*1.*JxW;
      (*ele_mat)(NDIM * i + 1, 16 ) += bas_grad[quad][i][2]*quad_coord[0]*JxW;
      (*ele_mat)(NDIM * i + 1, 17 ) += bas_grad[quad][i][2]*quad_coord[1]*JxW;
      (*ele_mat)(NDIM * i + 1, 18 ) += bas_grad[quad][i][2]*quad_coord[2]*JxW;
      (*ele_mat)(NDIM * i + 1, 19 ) += bas_grad[quad][i][2]*1.*JxW;
      /// 自由度i的z分量 653
      (*ele_mat)(NDIM * i + 2, 20 ) += bas_grad[quad][i][0]*quad_coord[0]*JxW;
      (*ele_mat)(NDIM * i + 2, 21 ) += bas_grad[quad][i][0]*quad_coord[1]*JxW;
      (*ele_mat)(NDIM * i + 2, 22 ) += bas_grad[quad][i][0]*quad_coord[2]*JxW;
      (*ele_mat)(NDIM * i + 2, 23 ) += bas_grad[quad][i][0]*1.*JxW;
      (*ele_mat)(NDIM * i + 2, 16 ) += bas_grad[quad][i][1]*quad_coord[0]*JxW;
      (*ele_mat)(NDIM * i + 2, 17 ) += bas_grad[quad][i][1]*quad_coord[1]*JxW;
      (*ele_mat)(NDIM * i + 2, 18 ) += bas_grad[quad][i][1]*quad_coord[2]*JxW;
      (*ele_mat)(NDIM * i + 2, 19 ) += bas_grad[quad][i][1]*1.*JxW;
      (*ele_mat)(NDIM * i + 2, 8 ) += bas_grad[quad][i][2]*quad_coord[0]*JxW;
      (*ele_mat)(NDIM * i + 2, 9 ) += bas_grad[quad][i][2]*quad_coord[1]*JxW;
      (*ele_mat)(NDIM * i + 2, 10 ) += bas_grad[quad][i][2]*quad_coord[2]*JxW;
      (*ele_mat)(NDIM * i + 2, 11 ) += bas_grad[quad][i][2]*1.*JxW;

    }

  }




}


/*************************************************************************
 * 计算单元刚度矩阵.
 *
 * update #3:材料相关部分
 ************************************************************************/
void LinearTet::buildStiffElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, double T_val) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  //update #3:
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");
  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");

  /*
    if(entity_id==1)
         material =	material_manager->getMaterial("Silicon");
    else if(entity_id==13)
         material =	material_manager->getMaterial("MoCu");
    else if(entity_id==2)
         material =	material_manager->getMaterial("Copper");
    else if(entity_id==4)
         material =	material_manager->getMaterial("SiO2");
    else if(entity_id==5)
         material =	material_manager->getMaterial("SiN");
    else if(entity_id==3)
         material =	material_manager->getMaterial("Gold");
    else if(entity_id==15)
         material =	material_manager->getMaterial("GaN");
    else if(entity_id==16)
         material =	material_manager->getMaterial("Al2O3");
    else if(entity_id==17)
         material =	material_manager->getMaterial("Alloy");
    else
         material =	material_manager->getMaterial("Gold");
                */

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_mat->resize(n_dof * NDIM, n_dof * NDIM);
  for (int i = 0; i < n_dof * NDIM; ++i) {
    for (int j = 0; j < n_dof * NDIM; ++j) {
      (*ele_mat)(i, j) = 0.0;
    }
  }
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);

  tbox::Array<tbox::Array<double> > moduli = material->getModuli(T_val);

  double a = moduli[0][0];
  double b = moduli[0][1];
  double c = moduli[3][3];
  /// 计算单元刚度矩阵.
  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      for (int l = 0; l < num_quad_pnts; ++l) {
        double JxW = volume * jac * weight[l];
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

/*************************************************************************
 * 计算单元质量矩阵.
 *
 * update #3:材料相关部分
 ************************************************************************/
void LinearTet::buildMassElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, double T_val) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  //update #3:
  //update #3:
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");
  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_mat->resize(n_dof * NDIM, n_dof * NDIM);
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

  double density =material->getDensity(T_val);

  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < NDIM; ++j) {
      (*ele_mat)(NDIM* i + j, NDIM * i + j) = 0.25 * volume * jac * density;
    }
  }
#if 0 
  ///  计算单元质量矩阵, 非集中质量.
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> >
      quad_pnt = integrator->getQuadraturePoints(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> >
      bas_val = shape_func->value(real_vertex,quad_pnt);
  for (int i = 0;i < n_dof; ++ i) {
    for (int j = 0;j < n_dof; ++ j) {
      for (int l = 0;l < num_quad_pnts;++ l) {
        double JxW = volume*jac*weight[l];
        for(int m = 0; m < NDIM; ++m){
          (*ele_mat)(NDIM*i+m,NDIM*j+m) +=
              JxW*density*(bas_val[l][i]*bas_val[l][j]);
        }
      }
    }
  }
#endif
}

/*************************************************************************
 * 计算单元阻尼矩阵.
 *
 * update #3:材料相关部分
 * //update #6 计算阻尼矩阵 at 2017-04-26
 ************************************************************************/
void LinearTet::buildDumpElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, double T_val)
{
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  //update #3:
  //update #3:
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");
  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_mat->resize(n_dof * NDIM, n_dof * NDIM);
  for(int i=0;i<n_dof*NDIM;i++)
    for(int j=0;j<n_dof*NDIM;j++)
      (*ele_mat)(i,j)=0;
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

  double damping =material->getDamping();
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> >
      quad_pnt = integrator->getQuadraturePoints(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> >
      bas_val = shape_func->value(real_vertex,quad_pnt);
  for (int i = 0;i < n_dof; ++ i) {
    for (int j = 0;j < n_dof; ++ j) {
      for (int l = 0;l < num_quad_pnts;++ l) {
        double JxW = volume*jac*weight[l];
        for(int m = 0; m < NDIM; ++m){
          (*ele_mat)(NDIM*i+m,NDIM*j+m) +=
              JxW*damping*(bas_val[l][i]*bas_val[l][j]);
        }
      }
    }
  }
}

/*************************************************************************
 * 计算单元矩阵.
 *
 * 计算单元矩阵P
 * P=M/(dt^2)+C/(2*dt)+beta*K
 *
 * //update #6  at 2017-05-02
 * Tested OK -05-05
 ************************************************************************/
void LinearTet::buildElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id,  tbox::Array<double> T_val)
{
  /// 取出积分器对象. 这里好像用不到
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();

  //计算单元刚度矩阵、质量矩阵和阻尼矩阵
  //刚度矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_K = new tbox::Matrix<double>();
  ele_mat_K->resize(n_dof* NDIM, n_dof* NDIM);
  //质量矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_M = new tbox::Matrix<double>();
  ele_mat_M->resize(n_dof* NDIM, n_dof* NDIM);
  //阻尼矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_C = new tbox::Matrix<double>();
  ele_mat_C->resize(n_dof* NDIM, n_dof* NDIM);

  //初始化
  for (int i1 = 0; i1 < n_dof* NDIM; ++i1) {
    for (int j = 0; j < n_dof* NDIM; ++j) {
      (*ele_mat_K)(i1, j) = 0.0;
      (*ele_mat_M)(i1, j) = 0.0;
      (*ele_mat_C)(i1, j) = 0.0;
    }
  }
  double e_Temperature=0;
  for(int i=0;i<4;i++)
    e_Temperature+=T_val[i]/4;
  //计算
  buildStiffElementMatrix(real_vertex, dt, time, ele_mat_K, entity_id,e_Temperature);
  buildMassElementMatrix(real_vertex, dt, time, ele_mat_M, entity_id,e_Temperature);
  buildDumpElementMatrix(real_vertex, dt, time, ele_mat_C, entity_id,e_Temperature);

  ele_mat->resize(n_dof * NDIM, n_dof * NDIM);
  for (int i = 0; i < n_dof * NDIM; ++i) {
    for (int j = 0; j < n_dof * NDIM; ++j) {
      (*ele_mat)(i, j) = 0.0;
    }
  }

  for (int i = 0; i < n_dof * NDIM; ++i) {
    for (int j = 0; j < n_dof * NDIM; ++j) {
      //	      (*ele_mat)(i, j) = (*ele_mat_M)(i, j)/(dt*dt)+(*ele_mat_C)(i, j)/(2*dt)+(*ele_mat_K)(i, j)*BETA;
      (*ele_mat)(i, j) = (*ele_mat_K)(i, j);
    }
  }
}

/*************************************************************************
 * 计算单元右端项..
 * Tested OK -05-05
 ************************************************************************/
void LinearTet::buildElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    NewmarkData *d_newmark,int entity_id,tbox::Array<double> T_val,
    tbox::Array<double> Tolder_val) {
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");
  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_vec->resize(n_dof * NDIM);
  for (int j = 0; j < n_dof * NDIM; ++j) {
    (*ele_vec)[j] = 0.0;
  }

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();
  //update #3:
  //update #3:
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");
  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");

  //update #8 计算单元温度 直接线性插值计算
  double e_Temperature=0;
  for(int i=0;i<n_dof;i++)
    e_Temperature+=T_val[i]/4;
  double e_Temperature_older=0;
  for(int i=0;i<n_dof;i++)
    e_Temperature_older+=Tolder_val[i]/4;

  //模量矩阵
  tbox::Array<tbox::Array<double> > moduli = material->getModuli(e_Temperature);
  double a = moduli[0][0];
  double b = moduli[0][1];
  //double c = moduli[3][3];

  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);

  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  //update #6 05-05pm
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

  //计算热应变
  //应变矩阵=th_stress[1 1 1 0 0 0]
  double alpha=material->getAlpha(e_Temperature);
  //  double th_stress=alpha*(e_Temperature-293.15);
  double th_stress=alpha*(e_Temperature-650.);
  //取出积分点处基函数梯度值
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);
  //热应变单元右端项
  double matrix_BTD[12];
  for(int j1=0;j1<12;j1++)
    matrix_BTD[j1]=0;
  //模量矩阵*热应力=termal_strain*[1 1 1 0 0 0]T
  double thermal_strain=(a+2*b)*th_stress;

  for(int i=0;i<n_dof;i++)
  {
    for(int j=0;j<NDIM;j++)
    {
      for(int l=0;l<num_quad_pnts;l++)
      {
        matrix_BTD[i*NDIM+j]+=thermal_strain*bas_grad[l][i][j]*weight[l]*jac*volume;
        //cout<<i*NDIM+j<<"   :"<<bas_grad[l][i][j]<<"   BTD:"<<matrix_BTD[i*NDIM+j]<<" thermal_strain*:"<<(thermal_strain*bas_grad[l][i][j]*weight[l]*jac)<<"  thermal_strain:"<<thermal_strain<<endl;
      }
      //  cout<<matrix_BTD[i*NDIM+j]<<"   ";
    }
    // cout<<'\n';
  }

  //静态时
#if 1

  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxWb = volume * jac * weight[l] * bas_val[l][i];
      for (int j = 0; j < NDIM; ++j) {
        (*ele_vec)[i * NDIM + j] += matrix_BTD[i*NDIM+j];
      }
    }
  }
#endif

  //瞬态时
#if 0
  //计算单元刚度矩阵、质量矩阵和阻尼矩阵
  //刚度矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_K = new tbox::Matrix<double>();
  ele_mat_K->resize(n_dof* NDIM, n_dof* NDIM);
  //质量矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_M = new tbox::Matrix<double>();
  ele_mat_M->resize(n_dof* NDIM, n_dof* NDIM);
  //阻尼矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_C = new tbox::Matrix<double>();
  ele_mat_C->resize(n_dof* NDIM, n_dof* NDIM);

  //初始化
  for (int i1 = 0; i1 < n_dof* NDIM; ++i1) {
    for (int j = 0; j < n_dof* NDIM; ++j) {
      (*ele_mat_K)(i1, j) = 0.0;
      (*ele_mat_M)(i1, j) = 0.0;
      (*ele_mat_C)(i1, j) = 0.0;
    }
  }

  //计算
  buildStiffElementMatrix(real_vertex, dt, time, ele_mat_K, entity_id,e_Temperature);
  buildMassElementMatrix(real_vertex, dt, time, ele_mat_M, entity_id,e_Temperature);
  buildDumpElementMatrix(real_vertex, dt, time, ele_mat_C, entity_id,e_Temperature);


  //计算Q、R矩阵
  //Q=-[-(2*M)/(Dt^2)+(1-2*beta)*K]
  //R=-[M/(Dt^2)-c/(2*Dt)+beta*K]
  //初始化矩阵Q、R
  tbox::Pointer<tbox::Matrix<double> > ele_mat_Q = new tbox::Matrix<double>();
  ele_mat_Q->resize(n_dof* NDIM, n_dof* NDIM);
  tbox::Pointer<tbox::Matrix<double> > ele_mat_R = new tbox::Matrix<double>();
  ele_mat_R->resize(n_dof* NDIM, n_dof* NDIM);

  for (int i1 = 0; i1 < n_dof* NDIM; ++i1) {
    for (int j = 0; j < n_dof* NDIM; ++j) {
      (*ele_mat_Q)(i1, j) = ((2*(*ele_mat_M)(i1, j))/(dt*dt))-((1-2*BETA)*(*ele_mat_K)(i1, j));
      (*ele_mat_R)(i1, j) = -(((*ele_mat_M)(i1, j)/(dt*dt))-((*ele_mat_C)(i1, j)/(2*dt))+(BETA*(*ele_mat_K)(i1, j)));
    }
  }

  //计算Q*solution（t-Dt）和R*solution(t-2*Dt)
  //定义一个向量用于存储计算结果、
  //update #6 -05-05  指针分配内存
  tbox::Pointer<tbox::Vector<double> > solutions_vec= new tbox::Vector<double>();
  solutions_vec->resize(n_dof * NDIM);
  for (int j = 0; j < n_dof * NDIM; ++j) {
    (*solutions_vec)[j] = 0.0;
  }

  for (int i1 = 0; i1 < n_dof* NDIM; ++i1) {
    for (int j = 0; j < n_dof* NDIM; ++j) {
      (*solutions_vec)[i1]+=(((*ele_mat_Q)(i1, j))*((d_newmark[j/3]).v_solution_old[j%3]))+
          (((*ele_mat_R)(i1, j))*((d_newmark[j/3]).v_solution_older[j%3]));
    }
  }

  //
  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxWb = volume * jac * weight[l] * bas_val[l][i];
      for (int j = 0; j < NDIM; ++j) {
        //   if (j % 3 == 2) {
        (*ele_vec)[i * NDIM + j] +=matrix_BTD[i*NDIM+j]-JxWb * 0.0e-3;
        //    }
      }
    }
  }
  for (int i1 = 0; i1 < n_dof; ++i1) {
    for(int j1=0;j1<NDIM;j1++)
      (*ele_vec)[i1*NDIM+j1]=BETA*((*ele_vec)[i1*NDIM+j1])+(1-2*BETA)*((d_newmark[i1]).v_rhs_old[j1])+
        BETA*((d_newmark[i1]).v_rhs_older[j1])+(*solutions_vec)[i1*NDIM+j1];
  }
#endif
}

/////////////////////////////////////////////////update #8//////////////////////////////////////////////

/*************************************************************************
 * 计算热计算单元矩阵.
 *
 ************************************************************************/
void LinearTet::buildTh_ElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, tbox::Array<double> T_val)
{
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  //update #3:
  //update #3:
  double e_Temperature=0;
  for(int i=0;i<4;i++)
    e_Temperature+=T_val[i]/4;
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");
  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");

  /// 取出单元上自由度数目.
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
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  double density=material->getDensity(e_Temperature);
  double K=material->getK(e_Temperature);
  double Cp=material->getCp(e_Temperature);

  /// 计算单元刚度矩阵.
  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      for (int l = 0; l < num_quad_pnts; ++l) {
        double JxW = volume * jac * weight[l];
        (*ele_mat)(i,j) += JxW*(density*Cp*bas_val[l][i]*bas_val[l][j]+
                                (bas_grad[l][i][0]*bas_grad[l][j][0]+
            bas_grad[l][i][1]*bas_grad[l][j][1]+
            bas_grad[l][i][2]*bas_grad[l][j][2])*K*dt);
      }
    }
  }
}


/*************************************************************************
 * 计算热计算单元矩阵.
 *
 ************************************************************************/
void LinearTet::buildStaticTh_ElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, tbox::Array<double> T_val)
{
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  //update #3:
  //update #3:
  double e_Temperature=0;
  for(int i=0;i<4;i++)
    e_Temperature+=T_val[i]/4;
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");
  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");

  /// 取出单元上自由度数目.
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
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  double K=material->getK(e_Temperature);

  /// 计算单元刚度矩阵.
  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      for (int l = 0; l < num_quad_pnts; ++l) {
        double JxW = volume * jac * weight[l];
        (*ele_mat)(i,j) += JxW*((bas_grad[l][i][0]*bas_grad[l][j][0]+
            bas_grad[l][i][1]*bas_grad[l][j][1]+
            bas_grad[l][i][2]*bas_grad[l][j][2])*K);
      }
    }
  }
}

/*************************************************************************
  * 计算热计算单元右端项.
  *
  ************************************************************************/
void LinearTet::buildTh_ElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    int entity_id, tbox::Array<double> T_val, double e_ThermalSource)
{
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();
  //update #3:
  //update #3:
  double e_Temperature=0;
  for(int i=0;i<4;i++)
    e_Temperature+=T_val[i]/4;
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");

  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_vec->resize(n_dof);
  for (int j = 0; j < n_dof; ++j) {
    (*ele_vec)[j] = 0.0;
  }


  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);

  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  //静态时
#if STATIC_S
  //update #6 05-05pm
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

  double f1=0.0;//热源
  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxWb = volume * jac * weight[l] * bas_val[l][i];
      (*ele_vec)[i] += JxWb * f1 * dt;
    }
  }
}
#endif

//瞬态时
#if TIME_S
//rhs=
//update #6 05-05pm
/// 取出积分点数目.
int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
/// 取出模板单元的面积.
double volume = integrator->getElementVolume();
/// 取出jacobian矩阵行列式.
double jac = integrator->getLocal2GlobalJacobian(real_vertex);

//先计算单元质量矩阵
tbox::Pointer<tbox::Matrix<double> > ele_mat_M = new tbox::Matrix<double>();
ele_mat_M->resize(n_dof,n_dof);
double density=material->getDensity(e_Temperature);
double Cp=material->getCp(e_Temperature);
for (int i = 0; i < n_dof; ++i) {
  for (int j = 0; j < n_dof; ++j) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxW = volume * jac * weight[l];
      (*ele_mat_M)(i,j) += JxW*density*Cp*(bas_val[l][i]*bas_val[l][j]);
    }
  }
}

/// M/dt T^{n-1}
for (int i1 = 0; i1 < n_dof; ++i1) {
  for (int j = 0; j < n_dof; ++j) {
    (*ele_vec)[i1]+=((*ele_mat_M)(i1, j))*T_val[j];
  }
}

/// 计算单元右端项.
for (int i = 0; i < n_dof; ++i) {
  for (int l = 0; l < num_quad_pnts; ++l) {
    double JxW = volume * jac * weight[l];
    (*ele_vec)[i] += JxW * e_ThermalSource * bas_val[l][i] * dt;
  }
}
#endif
}




/*************************************************************************
  * 计算热计算单元右端项(稳态).
  *
  ************************************************************************/
void LinearTet::buildStaticTh_ElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    int entity_id, tbox::Array<double> T_val, double e_ThermalSource)
{
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();
  //update #3:
  //update #3:
  double e_Temperature=0;
  for(int i=0;i<4;i++)
    e_Temperature+=T_val[i]/4;
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");

  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_vec->resize(n_dof);
  for (int j = 0; j < n_dof; ++j) {
    (*ele_vec)[j] = 0.0;
  }


  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);

  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);

  //瞬态时

  //rhs=
  //update #6 05-05pm
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

  //先计算单元质量矩阵
  tbox::Pointer<tbox::Matrix<double> > ele_mat_M = new tbox::Matrix<double>();
  ele_mat_M->resize(n_dof,n_dof);


  /// 计算单元右端项.
  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxW = volume * jac * weight[l];
      (*ele_vec)[i] += JxW * e_ThermalSource * bas_val[l][i];
    }
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////update #9//////////////////////////////////////////////

/*************************************************************************
  * 计算电计算单元矩阵.
  *
  ************************************************************************/
void LinearTet::buildE_ElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,int entity_id, tbox::Array<double> T_val)
{
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();

  //update #3:
  //update #3:
  double e_Temperature=0;
  for(int i=0;i<4;i++)
    e_Temperature+=T_val[i]/4;
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");

  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");
  /// 取出单元上自由度数目.
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
  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
      shape_func->gradient(real_vertex, quad_pnt);

  double sigma=material->getSigma(e_Temperature);

  /// 计算单元刚度矩阵.
  for (int i = 0; i < n_dof; ++i) {
    for (int j = 0; j < n_dof; ++j) {
      for (int l = 0; l < num_quad_pnts; ++l) {
        double JxW = volume * jac * weight[l];
        (*ele_mat)(i,j) += JxW*((bas_grad[l][i][0]*bas_grad[l][j][0]+
            bas_grad[l][i][1]*bas_grad[l][j][1]+
            bas_grad[l][i][2]*bas_grad[l][j][2])*sigma);
      }
    }
  }
}

/*************************************************************************
   * 计算电计算单元右端项.
   *
   ************************************************************************/
void LinearTet::buildE_ElementRHS(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Vector<double> > ele_vec,
    int entity_id, tbox::Array<double> T_val)
{
  /// 取出积分器对象.
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator =
      integrator_manager->getIntegrator("LinearTetrahedron");

  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");

  /// 取出材料.
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();
  //update #3:
  //update #3:
  double e_Temperature=0;
  for(int i=0;i<4;i++)
    e_Temperature+=T_val[i]/4;
  tbox::Pointer<Material> material =
      material_manager->getMaterial("Gold");
  if(entity_id==1)
    material =	material_manager->getMaterial("Silicon");
  else if(entity_id==13)
    material =	material_manager->getMaterial("MoCu");
  else if(entity_id==2)
    material =	material_manager->getMaterial("Copper");
  else if(entity_id==4)
    material =	material_manager->getMaterial("SiO2");
  else if(entity_id==5)
    material =	material_manager->getMaterial("SiN");
  else if(entity_id==3)
    material =	material_manager->getMaterial("Gold");
  else if(entity_id==15)
    material =	material_manager->getMaterial("GaN");
  else if(entity_id==16)
    material =	material_manager->getMaterial("Al2O3");
  else if(entity_id==17)
    material =	material_manager->getMaterial("Alloy");
  else if(entity_id==6)
    material =	material_manager->getMaterial("Aluminum");
  else
    material =	material_manager->getMaterial("Gold");

  /// 取出单元上自由度数目.
  int n_dof = shape_func->getNumberOfDof();
  ele_vec->resize(n_dof);
  for (int j = 0; j < n_dof; ++j) {
    (*ele_vec)[j] = 0.0;
  }


  /// 取出积分点.
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
      integrator->getQuadraturePoints(real_vertex);

  /// 取出积分点的积分权重.
  tbox::Array<double> weight = integrator->getQuadratureWeights();

  /// 取出基函数在积分点的值和梯度值.
  tbox::Array<tbox::Array<double> > bas_val =
      shape_func->value(real_vertex, quad_pnt);


  //update #6 05-05pm
  /// 取出积分点数目.
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  /// 取出模板单元的面积.
  double volume = integrator->getElementVolume();
  /// 取出jacobian矩阵行列式.
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

  double q=0.0;//电荷源
  for (int i = 0; i < n_dof; ++i) {
    for (int l = 0; l < num_quad_pnts; ++l) {
      double JxWb = volume * jac * weight[l] * bas_val[l][i];
      (*ele_vec)[i] += JxWb * q;
    }
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////

/*************************************************************************
 * 计算斯托克斯流 (Stokes Flow) 的初始单元刚度矩阵
 *
 * 物理背景: 用于流体非线性 Newton-Raphson 迭代的冷启动 (Cold Start)
 * 单元类型: P1+P1 (速度和压力都是线性插值)
 * 稳定化:   引入了 PSPG (压力稳定化彼得罗夫-伽辽金) 项以满足 LBB 条件
 * 矩阵结构: 节点交替自由度 [u1, v1, w1, p1, u2, v2, w2, p2 ...]
 ************************************************************************/
void LinearTet::buildInitFluidElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,
    int entity_id, tbox::Array<double> T_val,double *mu){
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
        IntegratorManager<NDIM>::getManager();
    tbox::Pointer<BaseIntegrator<NDIM> > integrator =
        integrator_manager->getIntegrator("LinearTetrahedron");

    tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
        ShapeFunctionManager<NDIM>::getManager();
    tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
        shape_manager->getShapeFunction("LinearTetrahedron");

    /// 2. 取出材料并计算平均温度
    tbox::Pointer<MaterialManager<NDIM> > material_manager =
        MaterialManager<NDIM>::getManager();
    double e_Temperature = 0;
    for(int i = 0; i < 4; i++) {
      e_Temperature += T_val[i] / 4.0;
    }
    tbox::Pointer<Material> material
        = material_manager->getMaterial(GET_USER_MAT(entity_id));
    /// 4. 单元与积分信息
      int n_nodes = shape_func->getNumberOfDof(); // 对于线性四面体，节点数为 4
      int n_dof_total = n_nodes * (NDIM + 1);     // 流体单元总自由度: 4 * (3 + 1) = 16

      ele_mat->resize(n_dof_total, n_dof_total);
      for (int i = 0; i < n_dof_total; ++i) {
        for (int j = 0; j < n_dof_total; ++j) {
          (*ele_mat)(i, j) = 0.0;
        }
      }
      double mu_value = (*mu);

      int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
      double volume = integrator->getElementVolume();
      double jac = integrator->getLocal2GlobalJacobian(real_vertex);

      tbox::Array<hier::DoubleVector<NDIM> > quad_pnt =
          integrator->getQuadraturePoints(real_vertex);
      tbox::Array<double> weight = integrator->getQuadratureWeights();
      tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad =
          shape_func->gradient(real_vertex, quad_pnt);
      tbox::Array<tbox::Array<double> > bas_val =
          shape_func->value(real_vertex, quad_pnt);
      /// 计算 PSPG 稳定化参数 tau
      /// 特征长度 h_e 粗略估计为体积的立方根
      double h_e = pow(volume, 1.0 / 3.0);
      double tau_pspg = (h_e * h_e) / (12.0 * mu_value); // 经典的 Stokes 稳定化时间尺度
      for (int i = 0; i < n_nodes; ++i) {
        for (int j = 0; j < n_nodes; ++j) {
          // 预先计算出局部块的行列索引，体现节点交替 [u, v, w, p]
          int row_u = i * (NDIM + 1) + 0;
          int row_v = i * (NDIM + 1) + 1;
          int row_w = i * (NDIM + 1) + 2;
          int row_p = i * (NDIM + 1) + 3;

          int col_u = j * (NDIM + 1) + 0;
          int col_v = j * (NDIM + 1) + 1;
          int col_w = j * (NDIM + 1) + 2;
          int col_p = j * (NDIM + 1) + 3;

          for (int l = 0; l < num_quad_pnts; ++l) {
            double JxW = volume * jac * weight[l];

            /** ----------------------------------------------------
            粘性扩散矩阵 K_mu: mu * \int (\nabla N_i \cdot \nabla N_j)
            ----------------------------------------------------**/
            double diff = JxW * mu_value * (bas_grad[l][i][0] * bas_grad[l][j][0] +
                bas_grad[l][i][1] * bas_grad[l][j][1] +
                bas_grad[l][i][2] * bas_grad[l][j][2]);

            // ----------------------------------------------------
            // B. 压力梯度矩阵 G: - \int (\nabla \cdot v_i) * P_j
            // ----------------------------------------------------
            double G_x = -JxW * bas_grad[l][i][0] * bas_val[l][j];
            double G_y = -JxW * bas_grad[l][i][1] * bas_val[l][j];
            double G_z = -JxW * bas_grad[l][i][2] * bas_val[l][j];

            // ----------------------------------------------------
            // C. 散度约束矩阵 D (乘以-1凑成对称鞍点系统: D = G^T)
            // 原本是 \int q_i * (\nabla \cdot u_j)，乘以-1后与 G_x 完全对称
            // ----------------------------------------------------
            double D_x = -JxW * bas_val[l][i] * bas_grad[l][j][0];
            double D_y = -JxW * bas_val[l][i] * bas_grad[l][j][1];
            double D_z = -JxW * bas_val[l][i] * bas_grad[l][j][2];

            // ----------------------------------------------------
            // D. PSPG 稳定化矩阵 K_pspg: \int tau * (\nabla q_i \cdot \nabla p_j)
            // 注意：因为连续性方程乘了-1，这里也必须乘-1，变成 -tau * ...
            // ----------------------------------------------------
            double pspg = -JxW * tau_pspg * (bas_grad[l][i][0] * bas_grad[l][j][0] +
                bas_grad[l][i][1] * bas_grad[l][j][1] +
                bas_grad[l][i][2] * bas_grad[l][j][2]);

            // ----------------------------------------------------
            // 组装到 16x16 的局部矩阵中 (NDIM=3 时)
            // ----------------------------------------------------

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
/*************************************************************************
 * 计算 N-S 方程的牛顿雅可比矩阵 (Newton Jacobian Matrix)
 *
 * 物理背景: 包含 Stokes 的纯粘性部分，外加对流项的严格偏导数 (偏导拆分)
 * 单元类型: P1+P1 (速度和压力都是线性插值), 包含 PSPG 稳定化
 * 额外输入: U_val - 从绘图量中提取的上一迭代步节点的已知速度向量 u^k
 ************************************************************************/
void LinearTet::buildFluidJacobianElementMatrix(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Pointer<tbox::Matrix<double> > ele_mat,
    int entity_id, tbox::Array<double> T_val,
    tbox::Array<hier::DoubleVector<NDIM> > U_val){
  /// 取出积分器与形函数对象
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager = IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator = integrator_manager->getIntegrator("LinearTetrahedron");
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager = ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func = shape_manager->getShapeFunction("LinearTetrahedron");
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();
  tbox::Pointer<Material> material
      = material_manager->getMaterial(GET_USER_MAT(entity_id));
  double e_Temperature = 0;
  for(int i = 0; i < 4; i++) {
    e_Temperature += T_val[i] / 4.0;
  }
  double mu = 1.002e-3;
  double rho = material->getDensity(e_Temperature);
  int n_nodes = shape_func->getNumberOfDof(); // 4
  int n_dof_total = n_nodes * (NDIM + 1);     // 16
  ele_mat->resize(n_dof_total, n_dof_total);
  for (int i = 0; i < n_dof_total; ++i) {
    for (int j = 0; j < n_dof_total; ++j) {
      (*ele_mat)(i, j) = 0.0;
    }
  }
  /// 提取积分点信息
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);

  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt = integrator->getQuadraturePoints(real_vertex);
  tbox::Array<double> weight = integrator->getQuadratureWeights();
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad = shape_func->gradient(real_vertex, quad_pnt);
  tbox::Array<tbox::Array<double> > bas_val = shape_func->value(real_vertex, quad_pnt);

  // 特征长度 h_e，用于稳定化参数计算
  double h_e = pow(volume, 1.0 / 3.0);
  for (int l = 0; l < num_quad_pnts; ++l){
    double JxW = volume * jac * weight[l];
    /// 计算当前积分点上的已知宏观速度 u^k 及其空间梯度 \nabla u^k
    double u_k = 0.0, v_k = 0.0, w_k = 0.0;
    double du_dx = 0.0, du_dy = 0.0, du_dz = 0.0;
    double dv_dx = 0.0, dv_dy = 0.0, dv_dz = 0.0;
    double dw_dx = 0.0, dw_dy = 0.0, dw_dz = 0.0;
    for (int m = 0; m < n_nodes; ++m) {
      double N_m = bas_val[l][m];
      double dNm_dx = bas_grad[l][m][0];
      double dNm_dy = bas_grad[l][m][1];
      double dNm_dz = bas_grad[l][m][2];

      u_k += N_m * U_val[m][0];
      v_k += N_m * U_val[m][1];
      w_k += N_m * U_val[m][2];

      du_dx += dNm_dx * U_val[m][0]; du_dy += dNm_dy * U_val[m][0]; du_dz += dNm_dz * U_val[m][0];
      dv_dx += dNm_dx * U_val[m][1]; dv_dy += dNm_dy * U_val[m][1]; dv_dz += dNm_dz * U_val[m][1];
      dw_dx += dNm_dx * U_val[m][2]; dw_dy += dNm_dy * U_val[m][2]; dw_dz += dNm_dz * U_val[m][2];
    }
    double U_norm = sqrt(u_k * u_k + v_k * v_k + w_k * w_k);

    /// PSPG 稳定化系数
    double tau_pspg = (h_e * h_e) / (12.0 * mu);

    /// SUPG 稳定化系数 (考虑对流与扩散的调和)
    double tau_supg = 0.0;
    if (U_norm > 1e-12) {
      tau_supg = 1.0 / sqrt(pow(2.0 * U_norm / h_e, 2.0) + pow(4.0 * mu / (rho * h_e * h_e), 2.0));
    } else {
      tau_supg = (h_e * h_e) / (4.0 * mu / rho);
    }
    /// 组装矩阵
    for(int i = 0; i < n_nodes; ++i){
      double Ni = bas_val[l][i];
      double dNi_dx = bas_grad[l][i][0];
      double dNi_dy = bas_grad[l][i][1];
      double dNi_dz = bas_grad[l][i][2];
      // 测试函数沿流线方向的偏导 (用于 SUPG)
      double U_dot_gradNi = u_k * dNi_dx + v_k * dNi_dy + w_k * dNi_dz;
      for (int j = 0; j < n_nodes; ++j){
        double Nj = bas_val[l][j];
        double dNj_dx = bas_grad[l][j][0];
        double dNj_dy = bas_grad[l][j][1];
        double dNj_dz = bas_grad[l][j][2];
        // 试探函数沿流线方向的偏导
        double U_dot_gradNj = u_k * dNj_dx + v_k * dNj_dy + w_k * dNj_dz;
        /// 单元自由度映射
        int row_u = i * 4 + 0; int col_u = j * 4 + 0;
        int row_v = i * 4 + 1; int col_v = j * 4 + 1;
        int row_w = i * 4 + 2; int col_w = j * 4 + 2;
        int row_p = i * 4 + 3; int col_p = j * 4 + 3;
        /// 粘性项
        double diff = JxW * mu * (dNi_dx * dNj_dx + dNi_dy * dNj_dy + dNi_dz * dNj_dz);

        double G_x = -JxW * dNi_dx * Nj;
        double G_y = -JxW * dNi_dy * Nj;
        double G_z = -JxW * dNi_dz * Nj;

        // 散度乘以 -1 (与 G 形成对称结构)
        double D_x = -JxW * Ni * dNj_dx;
        double D_y = -JxW * Ni * dNj_dy;
        double D_z = -JxW * Ni * dNj_dz;

        /// 标准 Galerkin 牛顿非线性对流雅可比
        /// 数学形式: \rho * Ni * (U^k \cdot \nabla \delta u) + \rho * Ni * (\delta u \cdot \nabla U^k)

        double rho_Ni = JxW * rho * Ni;

        // 主对角块 (U^k 传输 \delta u，外加 \delta u 自身方向的梯度扰动)
        double conv_uu = rho_Ni * (U_dot_gradNj + Nj * du_dx);
        double conv_vv = rho_Ni * (U_dot_gradNj + Nj * dv_dy);
        double conv_ww = rho_Ni * (U_dot_gradNj + Nj * dw_dz);

        // 非对角交叉块 (\delta u 在其他方向上的梯度扰动)
        double conv_uv = rho_Ni * Nj * du_dy;  // v_j 扰动 u 方程
        double conv_uw = rho_Ni * Nj * du_dz;  // w_j 扰动 u 方程

        double conv_vu = rho_Ni * Nj * dv_dx;  // u_j 扰动 v 方程
        double conv_vw = rho_Ni * Nj * dv_dz;  // w_j 扰动 v 方程

        double conv_wu = rho_Ni * Nj * dw_dx;  // u_j 扰动 w 方程
        double conv_wv = rho_Ni * Nj * dw_dy;  // v_j 扰动 w 方程
        /// ----------------------------------------------------
        /// SUPG & PSPG 稳定化项
        /// ----------------------------------------------------
        // PSPG 压力-压力项
        double pspg = -JxW * tau_pspg * (dNi_dx * dNj_dx + dNi_dy * dNj_dy + dNi_dz * dNj_dz);

        // SUPG 对流稳定项 (沿流线方向增加人工扩散)
        // \int \tau_{supg} (U^k \cdot \nabla N_i) * \rho * (U^k \cdot \nabla N_j)
        double supg_conv_base = JxW * tau_supg * rho * U_dot_gradNi;

        // SUPG 作用于牛顿雅可比的两部分偏导
        double supg_uu = supg_conv_base * (U_dot_gradNj + Nj * du_dx);
        double supg_vv = supg_conv_base * (U_dot_gradNj + Nj * dv_dy);
        double supg_ww = supg_conv_base * (U_dot_gradNj + Nj * dw_dz);

        double supg_uv = supg_conv_base * (Nj * du_dy);
        double supg_uw = supg_conv_base * (Nj * du_dz);
        double supg_vu = supg_conv_base * (Nj * dv_dx);
        double supg_vw = supg_conv_base * (Nj * dv_dz);
        double supg_wu = supg_conv_base * (Nj * dw_dx);
        double supg_wv = supg_conv_base * (Nj * dw_dy);

        // SUPG 压力梯度项
        // \int \tau_{supg} (U^k \cdot \nabla N_i) * \nabla N_j
        double supg_p_x = JxW * tau_supg * U_dot_gradNi * dNj_dx;
        double supg_p_y = JxW * tau_supg * U_dot_gradNi * dNj_dy;
        double supg_p_z = JxW * tau_supg * U_dot_gradNi * dNj_dz;


        /// ----------------------------------------------------
        /// 将所有项累加到 16x16 局部矩阵
        /// ----------------------------------------------------

        // --- 动量方程 X ---
        (*ele_mat)(row_u, col_u) += diff + conv_uu + supg_uu;
        (*ele_mat)(row_u, col_v) += conv_uv + supg_uv;
        (*ele_mat)(row_u, col_w) += conv_uw + supg_uw;
        (*ele_mat)(row_u, col_p) += G_x + supg_p_x;

        // --- 动量方程 Y ---
        (*ele_mat)(row_v, col_u) += conv_vu + supg_vu;
        (*ele_mat)(row_v, col_v) += diff + conv_vv + supg_vv;
        (*ele_mat)(row_v, col_w) += conv_vw + supg_vw;
        (*ele_mat)(row_v, col_p) += G_y + supg_p_y;

        // --- 动量方程 Z ---
        (*ele_mat)(row_w, col_u) += conv_wu + supg_wu;
        (*ele_mat)(row_w, col_v) += conv_wv + supg_wv;
        (*ele_mat)(row_w, col_w) += diff + conv_ww + supg_ww;
        (*ele_mat)(row_w, col_p) += G_z + supg_p_z;

        // --- 连续性方程 ---
        (*ele_mat)(row_p, col_u) += D_x;
        (*ele_mat)(row_p, col_v) += D_y;
        (*ele_mat)(row_p, col_w) += D_z;
        (*ele_mat)(row_p, col_p) += pspg;

      }
    }
  }

}
void LinearTet::buildFluidResidualElementVector(
    tbox::Array<hier::DoubleVector<NDIM> > real_vertex, const double dt,
    const double time, tbox::Array<double>& ele_vec,
    int entity_id, tbox::Array<hier::DoubleVector<NDIM> > U_val,
    tbox::Array<double> P_val,tbox::Array<double> T_val){
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager
      = IntegratorManager<NDIM>::getManager();
  tbox::Pointer<BaseIntegrator<NDIM> > integrator
      = integrator_manager->getIntegrator("LinearTetrahedron");
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager
      = ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func
      = shape_manager->getShapeFunction("LinearTetrahedron");
  tbox::Pointer<MaterialManager<NDIM> > material_manager =
      MaterialManager<NDIM>::getManager();
  tbox::Pointer<Material> material
      = material_manager->getMaterial(GET_USER_MAT(entity_id));
  double e_Temperature = 0;
  for(int i = 0; i < 4; i++) {
    e_Temperature += T_val[i] / 4.0;
  }
  double mu = 1.002e-3;
  double rho = material->getDensity(e_Temperature);

  int n_nodes = shape_func->getNumberOfDof(); // 4
  int num_quad_pnts = integrator->getNumberOfQuadraturePoints();
  double volume = integrator->getElementVolume();
  double jac = integrator->getLocal2GlobalJacobian(real_vertex);
  tbox::Array<hier::DoubleVector<NDIM> > quad_pnt = integrator
      ->getQuadraturePoints(real_vertex);
  tbox::Array<double> weight = integrator
      ->getQuadratureWeights();
  tbox::Array<tbox::Array<tbox::Array<double> > > bas_grad
      = shape_func->gradient(real_vertex, quad_pnt);
  tbox::Array<tbox::Array<double> > bas_val
      = shape_func->value(real_vertex, quad_pnt);

  double h_e = pow(volume, 1.0 / 3.0);
  for (int l = 0; l < num_quad_pnts; ++l){
    double JxW = volume * jac * weight[l];
    /// 积分点上的速度、压力和梯度变量
    double u_k = 0.0, v_k = 0.0, w_k = 0.0, p_k = 0.0;
    double du_dx = 0.0, du_dy = 0.0, du_dz = 0.0;
    double dv_dx = 0.0, dv_dy = 0.0, dv_dz = 0.0;
    double dw_dx = 0.0, dw_dy = 0.0, dw_dz = 0.0;
    double dp_dx = 0.0, dp_dy = 0.0, dp_dz = 0.0;
    for (int m = 0; m < n_nodes; ++m){
      double N_m = bas_val[l][m];
      double dNm_dx = bas_grad[l][m][0];
      double dNm_dy = bas_grad[l][m][1];
      double dNm_dz = bas_grad[l][m][2];

      u_k += N_m * U_val[m][0];
      v_k += N_m * U_val[m][1];
      w_k += N_m * U_val[m][2];
      p_k += N_m * P_val[m];

      du_dx += dNm_dx * U_val[m][0]; du_dy += dNm_dy * U_val[m][0]; du_dz += dNm_dz * U_val[m][0];
      dv_dx += dNm_dx * U_val[m][1]; dv_dy += dNm_dy * U_val[m][1]; dv_dz += dNm_dz * U_val[m][1];
      dw_dx += dNm_dx * U_val[m][2]; dw_dy += dNm_dy * U_val[m][2]; dw_dz += dNm_dz * U_val[m][2];

      dp_dx += dNm_dx * P_val[m]; dp_dy += dNm_dy * P_val[m]; dp_dz += dNm_dz * P_val[m];
    }
    /// =========================================================
    /// 计算强形式动量残差
    /// 注意: 对于线性四面体 (P1)，速度的二阶导数拉普拉斯项为 0，所以粘性残差消失。
    /// Mom_R = \rho (u \cdot \nabla u) + \nabla p
    /// =========================================================
    double Mom_Rx = rho * (u_k * du_dx + v_k * du_dy + w_k * du_dz) + dp_dx;
    double Mom_Ry = rho * (u_k * dv_dx + v_k * dv_dy + w_k * dv_dz) + dp_dy;
    double Mom_Rz = rho * (u_k * dw_dx + v_k * dw_dy + w_k * dw_dz) + dp_dz;

    /// 连续性强残差: div(u)
    double div_u = du_dx + dv_dy + dw_dz;

    /// =========================================================
    /// 步骤 C：计算稳定化参数
    /// =========================================================
    double U_norm = sqrt(u_k * u_k + v_k * v_k + w_k * w_k);
    double tau_pspg = (h_e * h_e) / (12.0 * mu);
    double tau_supg = 0.0;
    if (U_norm > 1e-12) {
      tau_supg = 1.0 / sqrt(pow(2.0 * U_norm / h_e, 2.0) + pow(4.0 * mu / (rho * h_e * h_e), 2.0));
    } else {
      tau_supg = (h_e * h_e) / (4.0 * mu / rho);
    }
    /// 组装积分
    for (int i = 0; i < n_nodes; ++i) {
      double Ni = bas_val[l][i];
      double dNi_dx = bas_grad[l][i][0];
      double dNi_dy = bas_grad[l][i][1];
      double dNi_dz = bas_grad[l][i][2];

      int row_u = i * 4 + 0;
      int row_v = i * 4 + 1;
      int row_w = i * 4 + 2;
      int row_p = i * 4 + 3;

      /// ----------------------------------------------------
      /// 1. 动量方程 (X, Y, Z) - 取【负残差】 (-R_u)
      /// ----------------------------------------------------
      // a. 纯粘性扩散积分: -\int \mu \nabla u : \nabla v
      double diff_u = -JxW * mu * (dNi_dx * du_dx + dNi_dy * du_dy + dNi_dz * du_dz);
      double diff_v = -JxW * mu * (dNi_dx * dv_dx + dNi_dy * dv_dy + dNi_dz * dv_dz);
      double diff_w = -JxW * mu * (dNi_dx * dw_dx + dNi_dy * dw_dy + dNi_dz * dw_dz);

      // b. 纯对流平流积分: -\int \rho (u \cdot \nabla u) v
      double conv_u = -JxW * Ni * rho * (u_k * du_dx + v_k * du_dy + w_k * du_dz);
      double conv_v = -JxW * Ni * rho * (u_k * dv_dx + v_k * dv_dy + w_k * dv_dz);
      double conv_w = -JxW * Ni * rho * (u_k * dw_dx + v_k * dw_dy + w_k * dw_dz);

      // c. 压力梯度积分: +\int p (\nabla \cdot v)  (弱形式分部积分翻转符号)
      double gradp_u = JxW * p_k * dNi_dx;
      double gradp_v = JxW * p_k * dNi_dy;
      double gradp_w = JxW * p_k * dNi_dz;

      // d. SUPG 稳定化负残差: -\int \tau_{supg} (u \cdot \nabla N_i) * (Mom_R)
      double U_dot_gradNi = u_k * dNi_dx + v_k * dNi_dy + w_k * dNi_dz;
      double supg_u = -JxW * tau_supg * U_dot_gradNi * Mom_Rx;
      double supg_v = -JxW * tau_supg * U_dot_gradNi * Mom_Ry;
      double supg_w = -JxW * tau_supg * U_dot_gradNi * Mom_Rz;

      /// 累加动量右端项
      ele_vec[row_u] += diff_u + conv_u + gradp_u + supg_u;
      ele_vec[row_v] += diff_v + conv_v + gradp_v + supg_v;
      ele_vec[row_w] += diff_w + conv_w + gradp_w + supg_w;

      /// ----------------------------------------------------
      /// 2. 连续性方程 (压力 P) - 取【正残差】 (+R_p)
      /// ----------------------------------------------------
      // a. 散度积分: +\int q (\nabla \cdot u)
      double div_term = JxW * Ni * div_u;

      // b. PSPG 稳定化正残差: +\int \tau_{pspg} (\nabla q \cdot Mom_R)
      double pspg_term = JxW * tau_pspg * (dNi_dx * Mom_Rx + dNi_dy * Mom_Ry + dNi_dz * Mom_Rz);

      /// 累加连续性方程右端项
      ele_vec[row_p] += div_term + pspg_term;
    }
  }

}

const int LinearTet::getProblemDim() { return NDIM; }

//获取实体自由度
tbox::Array<int> LinearTet::getNumberOfDofOnEntity() {
  /// 取出形函数对象.
  tbox::Pointer<ShapeFunctionManager<NDIM> > shape_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<BaseShapeFunction<NDIM> > shape_func =
      shape_manager->getShapeFunction("LinearTetrahedron");
  return shape_func->getNumberOfDofOnEntity();
}
