#include "CafeParsing.h"
#include "MapVoxelManager.h"
#include "FiniteElementManager.h"
#include "CellularAutomataManager.h"
#include "CellularAutomataManager_3D.h"
#include "Node.h"
#include "Element.h"
#include <yaml-cpp/yaml.h>
#include <assert.h>
#define __max(a,b)  (((a) > (b)) ? (a) : (b))
#define __min(a,b)  (((a) < (b)) ? (a) : (b))
Bin::Bin()
  : hasElements_(false)
{
  //nothing to do
}

Bin::~Bin()
{
  // nothing to do
}

/*----------------------
  ~ constructor
  ----------------------*/
MapVoxelManager::MapVoxelManager
(FiniteElementManager* feManager, CellularAutomataManager* caManager)
  : feManager_(feManager),
    caManager_(caManager),
    maxNonlinearIter_(20),
    isInElemTolerance_(1.0e-16),
    nullMatTemperature_(-1.0e8)
{
  int count = 0;
  for (int i = 0; i < 3; i++)
  {
    // Notation: tricky code count++
    xyzMinMax_[count++] =  1.0e8; // set min as this value
    xyzMinMax_[count++] = -1.0e8; // set max as this value

    sizeOfBin[i] = 10;
  }

}

/*------------------
    ~ destructor
  ------------------*/
MapVoxelManager::~MapVoxelManager()
{
  if (elementTagList_)
    delete[] elementTagList_;

  if (parCoordsList_)
    delete[] parCoordsList_;

}

/*---------------
   ~ load
  ---------------*/
void
MapVoxelManager::load(const YAML::Node& node)
{
  const YAML::Node * mapVoxel = node.FindValue("transfers");
  CafeEnv::self().caOutputP0() << "\n" << "Transfer Review Review" << "\n";
  CafeEnv::self().caOutputP0() << "=============================" << std::endl;

  if (mapVoxel)
  {
    for (size_t iMapVoxel = 0; iMapVoxel < mapVoxel->size(); iMapVoxel++)
    {
      const YAML::Node & mV_node = (*mapVoxel)[iMapVoxel];
      std::string name = "mapping_voxel";
      get_if_present_no_default(mV_node, "bin_size", sizeOfBin);
      get_if_present(mV_node, "tolerance", isInElemTolerance_, isInElemTolerance_);
      get_if_present(mV_node, "max_iterations", maxNonlinearIter_, maxNonlinearIter_);
      get_if_present(mV_node, "void_temperature", nullMatTemperature_, nullMatTemperature_);
      feManager_->nullMatTemperature_ = nullMatTemperature_;

      // yaml-cpp Receipt
      CafeEnv::self().caOutputP0() << name << "\n";
      CafeEnv::self().caOutputP0() << "Bin Size: " << sizeOfBin[0] << "\n";
      CafeEnv::self().caOutputP0() << "N-R Method Tolerance: " << isInElemTolerance_ << "\n";
      CafeEnv::self().caOutputP0() << "N-R Max-iteration: " << maxNonlinearIter_ << "\n";
      CafeEnv::self().caOutputP0() << std::endl;;
    }
  }
}

/*-------------------
   ~ initialize
  -------------------*/
void
MapVoxelManager::initialize()
{
  nDim_ = caManager_->nDim_;
  numVoxels_ = caManager_->numLocGhostCells_;
  elementTagList_ = new int[numVoxels_];
  parCoordsList_ = new double[nDim_*numVoxels_];

  for (int i = 0; i < numVoxels_; i++)
  {
    // assign minus value to indicate no element tag
    elementTagList_[i] = -1;   
    for (int j=0; j<nDim_; j++)
      parCoordsList_[i*nDim_ + j] = -2;
  }
}

/*----------------------------
   ~ execute
  ----------------------------*/
void
MapVoxelManager::execute()
{  
  // Step 1: get domain size coverd by fintite element mesh used in feManager
  get_domain_xyz_min_max();

  // Step 2: list of bins to improve searching efficiency
  std::vector<Bin> binList;
  int numBin[3];
  double dBin[3];
  double dVoxel = caManager_->h_;

  for (int i = 0; i < 3; i++)
  {
    dBin[i] = sizeOfBin[i] * dVoxel;
    numBin[i] = ceil((xyzMinMax_[i * 2 + 1] - xyzMinMax_[i * 2]) / dBin[i]);
  }
  binList.resize(numBin[0] * numBin[1] * numBin[2]);

  // Step 3: attache the elements to the bins
  attach_elements(dBin, numBin, binList);

  // Step 4: get the element-tag for voxles
  tag_voxels(dBin, numBin, binList);
}

void
MapVoxelManager::execute_with_focus_on_CA_region()
{
  // Step 1: get domain size coverd by fintite element mesh used in feManager
  get_CA_domain_xyz_min_max();

  // Step 2: list of bins to improve searching efficiency
  std::vector<Bin> binList;
  int numBin[3];
  double dBin[3];
  double dVoxel = caManager_->h_;

  for (int i = 0; i < 3; i++)
  {
    dBin[i] = sizeOfBin[i] * dVoxel;
    numBin[i] = ceil((xyzMinMax_[i * 2 + 1] - xyzMinMax_[i * 2]) / dBin[i]);
  }
  binList.resize(numBin[0] * numBin[1] * numBin[2]);

  // Step 3: attache the elements to the bins
  attach_elements(dBin, numBin, binList);

  // Step 4: get the element-tag for voxles
  tag_voxels(dBin, numBin, binList);
}
/*-----------------------------
    ~ update
  -----------------------------*/
void 
MapVoxelManager::update()
{
  // nothing for now
}

/*----------------------------
   find the max/mins of mesh
   ~ get_domain_xyz_min_max
  ----------------------------*/
void
MapVoxelManager::get_domain_xyz_min_max()
{
  int numbNodes = feManager_->numNodes_;
  for (int iNode = 0; iNode < numbNodes; iNode++)
  {
    double* coord = feManager_->nodeVector_[iNode].coordinates_;
    
    int count = 0;
    for (int iCom = 0; iCom < nDim_; iCom++)
    {
      // Notation: tricky code xyzMinMax_[count++] does not work on linux
      xyzMinMax_[count] = __min(xyzMinMax_[count], coord[iCom]);
      count++;
      xyzMinMax_[count] = __max(xyzMinMax_[count], coord[iCom]);
      count++;
    }
  }
}

void
MapVoxelManager::get_CA_domain_xyz_min_max()
{
  // Step 1: get the domain size of the CA region of the current process
  double dVoxel = caManager_->h_;  
  xyzMinMax_[0] = (caManager_->xStart_ - 1) * dVoxel + caManager_->x0_;
  xyzMinMax_[2] = (caManager_->yStart_ - 1) * dVoxel + caManager_->y0_;
  xyzMinMax_[4] = (caManager_->zStart_ - 1) * dVoxel + caManager_->z0_;

  xyzMinMax_[1] = xyzMinMax_[0] + caManager_->nxLocGhost_*dVoxel;
  xyzMinMax_[3] = xyzMinMax_[2] + caManager_->nyLocGhost_*dVoxel;
  xyzMinMax_[5] = xyzMinMax_[4] + caManager_->nzLocGhost_*dVoxel;

  // Step 2: initial the record of the current domain size
  double xyzMinMaxCurrent[6];
  for (int iCom = 0; iCom < 6; iCom++) {
    xyzMinMaxCurrent[iCom] = xyzMinMax_[iCom];
  }

  // Step 3: loop over all elements to determine the maximum domain size
  int numbElements = feManager_->numElements_; // global numbers
  for (int iEle = 0; iEle < numbElements; iEle++) 
  {
    Element* ele = feManager_->elementVector_[iEle];
    bool overlap = false;
    int numbNode = ele->get_nodes_per_element();
    for (int iNode = 0; iNode < numbNode; iNode++) {
      int nodeID = ele->nID_[iNode];
      double* coord = feManager_->nodeVector_[nodeID].coordinates_;
      if (coord[0] >= xyzMinMax_[0] && coord[0] <= xyzMinMax_[1] &&
        coord[1] >= xyzMinMax_[2] && coord[1] <= xyzMinMax_[3] &&
        coord[2] >= xyzMinMax_[4] && coord[2] <= xyzMinMax_[5]) {
        overlap = true;
        break;
      }
    }

    // expand the domain size
    if (overlap) 
    {
      double range[6];
      get_cell_xyz_min_max(range, feManager_->elementVector_[iEle]);
      
      int count = 0;
      for (int iCom = 0; iCom < nDim_; iCom++) {
        xyzMinMaxCurrent[count] = __min(xyzMinMaxCurrent[count], range[count]);
        count++;
        xyzMinMaxCurrent[count] = __max(xyzMinMaxCurrent[count], range[count]);
        count++;
      }
    }
  }

  // Step 4: reset the xyzMinMax_ as the final xyzMinMaxCurrent
  int count = 0;
  for (int iCom = 0; iCom < 6; iCom++) {
    xyzMinMax_[iCom] = xyzMinMaxCurrent[iCom];    
  }
}

/*----------------------------
   find the max/mins of mesh
  ~ get_cell_xyz_min_max
----------------------------*/
void
MapVoxelManager::get_cell_xyz_min_max(double * range, Element* ele)
{
  int count = 0;
  for (int i = 0; i < 3; i++)
  {
    // Notation: tricky code count++
    range[count++] =  1.0e8; // set min as this value
    range[count++] = -1.0e8; // set max as this value
  }

  int numNode = ele->get_nodes_per_element();
  for (int iNode = 0; iNode < numNode; iNode++)
  {
    int ID = ele->nID_[iNode];
    double* coord = feManager_->nodeVector_[ID].coordinates_;

    int count = 0;
    for (int iCom = 0; iCom < nDim_; iCom++)
    {
      // Notation: tricky code range[count++] does not work on linux
      range[count] = __min(range[count], coord[iCom]);
      count++;
      range[count] = __max(range[count], coord[iCom]);
      count++;
    }
  }
}

/*------------------------------------------
     attach elements to corresponding bins
  ------------------------------------------*/
void
MapVoxelManager::attach_elements(double* dBin, int* numBin, std::vector<Bin> &binList)
{
  int numElement = feManager_->numElements_;
  int numBinX = numBin[0];
  int numBinXY = numBin[0]* numBin[1];

  for (int iEle = 0; iEle < numElement; iEle++)
  {
    int elementID = iEle;  // FIXME: may use a pointer array to store the element ID
    double range[6];
    Element * ele = feManager_->elementVector_[iEle];
    get_cell_xyz_min_max(range, feManager_->elementVector_[iEle]);
    // the element should be wholely inside the domain defined by xyzMinMax
    bool inside = true;
    for (int iCom = 0; iCom < 3; iCom++) 
    {
      if (range[iCom * 2] < xyzMinMax_[iCom * 2] ||
          range[iCom * 2 + 1] > xyzMinMax_[iCom * 2 + 1]) {
        inside = false;
        break;
      }      
    }

    if (!inside) {
      continue;
    }

    // set nodal value
    int numNode = ele->get_nodes_per_element();
    for (int iNode = 0; iNode < numNode; iNode++)
    {
      int ID = ele->nID_[iNode];
      feManager_->nodeVector_[ID].grainAround_ = true;
    }

    // find bounds for bins
    int imin = std::floor((range[0] - xyzMinMax_[0]) / dBin[0]);
    int imax = std::ceil ((range[1] - xyzMinMax_[0]) / dBin[0]);
    int jmin = std::floor((range[2] - xyzMinMax_[2]) / dBin[1]);
    int jmax = std::ceil ((range[3] - xyzMinMax_[2]) / dBin[1]);
    int kmin = std::floor((range[4] - xyzMinMax_[4]) / dBin[2]);
    int kmax = std::ceil ((range[5] - xyzMinMax_[4]) / dBin[2]);

    for (int kk = kmin; kk < kmax; kk++)
    {
      for (int jj = jmin; jj < jmax; jj++)
      {
        for (int ii = imin; ii < imax; ii++)
        {
          int ind = ii + jj * numBinX + kk * numBinXY;
          binList[ind].hasElements_ = true;
          binList[ind].elements_.push_back(elementID);
        }//end for(ii)
      }//end for(jj)
    }//end for(kk)
  }

}

/*-------------------------------
    tag_voxels
  -------------------------------*/
void
MapVoxelManager::tag_voxels(double* dBin, int* numBin, std::vector<Bin> &binList)
{
  // prepare parameters for calculating the coordinates of the voxel
  double dVoxel = caManager_->h_;
  double originCoords[3];
  originCoords[0] = (caManager_->xStart_ - 0.5) * dVoxel + caManager_->x0_;
  originCoords[1] = (caManager_->yStart_ - 0.5) * dVoxel + caManager_->y0_;
  originCoords[2] = (caManager_->zStart_ - 0.5) * dVoxel + caManager_->z0_;
  bool findEle = false;
  int nxVoxel = caManager_->nxLocGhost_;
  int nxyVoxel = caManager_->nxyLocGhost_;

  // assign voxels to elements
  for (int iVoxel = 0; iVoxel < numVoxels_; iVoxel++)
  {
    // get the coordinates of iVoxel the the index of bins
    int indexVoxel[3];  // I J K;
    int indexMinMaxBins[6]; // Min-Max

    indexVoxel[2] = iVoxel / nxyVoxel;
    indexVoxel[1] = (iVoxel % nxyVoxel) / nxVoxel;
    indexVoxel[0] = (iVoxel % nxyVoxel) % nxVoxel;
    double coordVoxel[3];
    bool outOfRegion = false;

    for (int iCom = 0; iCom < 3; iCom++)
    {
      coordVoxel[iCom] = originCoords[iCom] + indexVoxel[iCom] * dVoxel;
      int min = 2 * iCom;
      int max = min + 1;
      if (coordVoxel[iCom] < xyzMinMax_[min] || coordVoxel[iCom] > xyzMinMax_[max])
      {
        outOfRegion = true;
        break;
      }
      double distance = coordVoxel[iCom] - xyzMinMax_[min];
      indexMinMaxBins[min] = std::floor(distance / dBin[iCom]);
      indexMinMaxBins[max] = std::ceil (distance / dBin[iCom]);
    }
    // out of region covered by finite element
    if (outOfRegion)
      continue;

    // adjust the idnex if necessary
    for (int i = 0; i < 3; i++)
    {
      int min = 2 * i;
      int max = min + 1;
      if (indexMinMaxBins[min] == indexMinMaxBins[max] &&
        indexMinMaxBins[min] > 0 &&
        indexMinMaxBins[max] < numBin[i])
      {
        indexMinMaxBins[min] -= 1;
        indexMinMaxBins[max] += 1;
      }
      else if (indexMinMaxBins[max] == 0)
      {
        indexMinMaxBins[max] = 1;
      }
      else if (indexMinMaxBins[max] == numBin[i])
      {
        if (indexMinMaxBins[min] > 0)
          indexMinMaxBins[min] -= 1;
      }
    }

    // tag element to iVoxel
    int imin = indexMinMaxBins[0];
    int imax = indexMinMaxBins[1];
    int jmin = indexMinMaxBins[2];
    int jmax = indexMinMaxBins[3];
    int kmin = indexMinMaxBins[4];
    int kmax = indexMinMaxBins[5];
    int binNx = numBin[0];
    int binNxy = numBin[0] * numBin[1];
    bool gotTheElement = false;
    for (int kk = kmin; kk < kmax; kk++)
    {
      for (int jj = jmin; jj < jmax; jj++)
      {
        for (int ii = imin; ii < imax; ii++)
        {
          int ind = ii + jj * binNx + kk * binNxy;

          if (binList[ind].hasElements_)
          {
            int numEle = binList[ind].elements_.size();
            std::vector<int>* binElements = &binList[ind].elements_;
            for (int ll = 0; ll < numEle; ll++)
            {
              int elementID = (*binElements)[ll];
              Element * ele = feManager_->elementVector_[elementID];
              double range[6];
              get_cell_xyz_min_max(range, ele);
              if (coordVoxel[0] < range[0] || coordVoxel[0] > range[1] ||
                  coordVoxel[1] < range[2] || coordVoxel[1] > range[3] ||
                  coordVoxel[2] < range[4] || coordVoxel[2] > range[5])
                continue;

              double paraCoords[3];              
              if (is_in_element(ele, coordVoxel, paraCoords))
              {
                elementTagList_[iVoxel] = elementID;
                for (int mm = 0; mm < 3; mm++)
                  parCoordsList_[iVoxel * 3 + mm] = paraCoords[mm];
                gotTheElement = true;
                findEle = true;
              }
              if (gotTheElement)
                break;
            } // end ll
          } // end current bin
          if (gotTheElement)
            break;
        } // end ii
        if (gotTheElement)
          break;
      } // end jj
      if (gotTheElement)
        break;
    } // end kk
  } // end iVoxel
  if (!findEle)
      throw std::runtime_error("Error: a proc has no corrsponding FEM mesh!! Maybe decrease the number of multi-core parallel can help.");
}

/*------------------------------
   ~ is_in_element
  ------------------------------*/
bool
MapVoxelManager::is_in_element(Element* ele, double *pointCoords, double* paraCoords)
{
  // Translate element so that (x,y,z) coordinates of the first node are (0,0,0)
  int numNode = ele->get_nodes_per_element();
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  x.resize(numNode);
  y.resize(numNode);
  z.resize(numNode);
  //double *x = new double[numNode];
  //double *y = new double[numNode];
  //double *z = new double[numNode];   
  x[0] = y[0] = z[0] = 0.0;

  int ID = ele->nID_[0];
  double *originPoint = feManager_->nodeVector_[ID].coordinates_;
  
  for (int iNode = 1; iNode < numNode; iNode++)
  {
    int ID = ele->nID_[iNode];
    double* coords = feManager_->nodeVector_[ID].coordinates_;
    x[iNode] = 0.125*(coords[0] - originPoint[0]);
    y[iNode] = 0.125*(coords[1] - originPoint[1]);
    z[iNode] = 0.125*(coords[2] - originPoint[2]);
  }

  // (xp,yp,zp) is the point at which we're searching for (xi,eta,zeta)
  // (must translate this also)

  double xp = pointCoords[0] - originPoint[0];
  double yp = pointCoords[1] - originPoint[1];
  double zp = pointCoords[2] - originPoint[2];

  // Newton-Raphson iteration for (xi,eta,zeta)
  double j[9];
  double f[3];
  double shapefct[8];
  double xinew = 0.5;     // initial guess
  double etanew = 0.5;
  double zetanew = 0.5;
  double xicur = 0.5;
  double etacur = 0.5;
  double zetacur = 0.5;
  double xidiff[] = { 1.0, 1.0, 1.0 };
  
  int i = 0;

  do
  {
    j[0] =
      -(1.0 - etacur)*(1.0 - zetacur)*x[1]
      - (1.0 + etacur)*(1.0 - zetacur)*x[2]
      + (1.0 + etacur)*(1.0 - zetacur)*x[3]
      + (1.0 - etacur)*(1.0 + zetacur)*x[4]
      - (1.0 - etacur)*(1.0 + zetacur)*x[5]
      - (1.0 + etacur)*(1.0 + zetacur)*x[6]
      + (1.0 + etacur)*(1.0 + zetacur)*x[7];

    j[1] =
      (1.0 + xicur)*(1.0 - zetacur)*x[1]
      - (1.0 + xicur)*(1.0 - zetacur)*x[2]
      - (1.0 - xicur)*(1.0 - zetacur)*x[3]
      + (1.0 - xicur)*(1.0 + zetacur)*x[4]
      + (1.0 + xicur)*(1.0 + zetacur)*x[5]
      - (1.0 + xicur)*(1.0 + zetacur)*x[6]
      - (1.0 - xicur)*(1.0 + zetacur)*x[7];

    j[2] =
        (1.0 - etacur)*(1.0 + xicur)*x[1]
      + (1.0 + etacur)*(1.0 + xicur)*x[2]
      + (1.0 + etacur)*(1.0 - xicur)*x[3]
      - (1.0 - etacur)*(1.0 - xicur)*x[4]
      - (1.0 - etacur)*(1.0 + xicur)*x[5]
      - (1.0 + etacur)*(1.0 + xicur)*x[6]
      - (1.0 + etacur)*(1.0 - xicur)*x[7];

    j[3] =
      - (1.0 - etacur)*(1.0 - zetacur)*y[1]
      - (1.0 + etacur)*(1.0 - zetacur)*y[2]
      + (1.0 + etacur)*(1.0 - zetacur)*y[3]
      + (1.0 - etacur)*(1.0 + zetacur)*y[4]
      - (1.0 - etacur)*(1.0 + zetacur)*y[5]
      - (1.0 + etacur)*(1.0 + zetacur)*y[6]
      + (1.0 + etacur)*(1.0 + zetacur)*y[7];

    j[4] =
        (1.0 + xicur)*(1.0 - zetacur)*y[1]
      - (1.0 + xicur)*(1.0 - zetacur)*y[2]
      - (1.0 - xicur)*(1.0 - zetacur)*y[3]
      + (1.0 - xicur)*(1.0 + zetacur)*y[4]
      + (1.0 + xicur)*(1.0 + zetacur)*y[5]
      - (1.0 + xicur)*(1.0 + zetacur)*y[6]
      - (1.0 - xicur)*(1.0 + zetacur)*y[7];

    j[5] =
        (1.0 - etacur)*(1.0 + xicur)*y[1]
      + (1.0 + etacur)*(1.0 + xicur)*y[2]
      + (1.0 + etacur)*(1.0 - xicur)*y[3]
      - (1.0 - etacur)*(1.0 - xicur)*y[4]
      - (1.0 - etacur)*(1.0 + xicur)*y[5]
      - (1.0 + etacur)*(1.0 + xicur)*y[6]
      - (1.0 + etacur)*(1.0 - xicur)*y[7];

    j[6] =
      - (1.0 - etacur)*(1.0 - zetacur)*z[1]
      - (1.0 + etacur)*(1.0 - zetacur)*z[2]
      + (1.0 + etacur)*(1.0 - zetacur)*z[3]
      + (1.0 - etacur)*(1.0 + zetacur)*z[4]
      - (1.0 - etacur)*(1.0 + zetacur)*z[5]
      - (1.0 + etacur)*(1.0 + zetacur)*z[6]
      + (1.0 + etacur)*(1.0 + zetacur)*z[7];

    j[7] =
        (1.0 + xicur)*(1.0 - zetacur)*z[1]
      - (1.0 + xicur)*(1.0 - zetacur)*z[2]
      - (1.0 - xicur)*(1.0 - zetacur)*z[3]
      + (1.0 - xicur)*(1.0 + zetacur)*z[4]
      + (1.0 + xicur)*(1.0 + zetacur)*z[5]
      - (1.0 + xicur)*(1.0 + zetacur)*z[6]
      - (1.0 - xicur)*(1.0 + zetacur)*z[7];

    j[8] =
        (1.0 - etacur)*(1.0 + xicur)*z[1]
      + (1.0 + etacur)*(1.0 + xicur)*z[2]
      + (1.0 + etacur)*(1.0 - xicur)*z[3]
      - (1.0 - etacur)*(1.0 - xicur)*z[4]
      - (1.0 - etacur)*(1.0 + xicur)*z[5]
      - (1.0 + etacur)*(1.0 + xicur)*z[6]
      - (1.0 + etacur)*(1.0 - xicur)*z[7];

    double jdet = -(j[2] * j[4] * j[6]) + j[1] * j[5] * j[6] + j[2] * j[3] * j[7] 
                  -j[0] * j[5] * j[7] - j[1] * j[3] * j[8] + j[0] * j[4] * j[8];

    if (!jdet) {
      i = maxNonlinearIter_;
      break;
    }
    shapefct[0] = (1.0 - etacur)*(1.0 - xicur)*(1.0 - zetacur);

    shapefct[1] = (1.0 - etacur)*(1.0 + xicur)*(1.0 - zetacur);

    shapefct[2] = (1.0 + etacur)*(1.0 + xicur)*(1.0 - zetacur);

    shapefct[3] = (1.0 + etacur)*(1.0 - xicur)*(1.0 - zetacur);

    shapefct[4] = (1.0 - etacur)*(1.0 - xicur)*(1.0 + zetacur);

    shapefct[5] = (1.0 - etacur)*(1.0 + xicur)*(1.0 + zetacur);

    shapefct[6] = (1.0 + etacur)*(1.0 + xicur)*(1.0 + zetacur);

    shapefct[7] = (1.0 + etacur)*(1.0 - xicur)*(1.0 + zetacur);

    f[0] = xp - shapefct[1] * x[1] - shapefct[2] * x[2] - shapefct[3] * x[3] - shapefct[4] * x[4] - 
                shapefct[5] * x[5] - shapefct[6] * x[6] - shapefct[7] * x[7];

    f[1] = yp - shapefct[1] * y[1] - shapefct[2] * y[2] - shapefct[3] * y[3] - shapefct[4] * y[4] - 
                shapefct[5] * y[5] - shapefct[6] * y[6] - shapefct[7] * y[7];

    f[2] = zp - shapefct[1] * z[1] - shapefct[2] * z[2] - shapefct[3] * z[3] - shapefct[4] * z[4] - 
                shapefct[5] * z[5] - shapefct[6] * z[6] - shapefct[7] * z[7];

    xinew = (jdet*xicur + f[2] * (j[2] * j[4] - j[1] * j[5]) 
                        + f[1] * (j[1] * j[8] - j[2] * j[7])
                        + f[0] * (j[5] * j[7] - j[4] * j[8]) 
                        )/jdet;

    etanew = (etacur*jdet + f[2] * (j[0] * j[5] - j[2] * j[3])
                          + f[1] * (j[2] * j[6] - j[0] * j[8])
                          + f[0] * (j[3] * j[8] - j[5] * j[6])
                          )/jdet;

    zetanew = (jdet*zetacur + f[2] * (j[1] * j[3] - j[0] * j[4]) 
                            + f[1] * (j[0] * j[7] - j[1] * j[6])
                            + f[0] * (j[4] * j[6] - j[3] * j[7])
                            )/jdet;

    xidiff[0] = xinew - xicur;
    xidiff[1] = etanew - etacur;
    xidiff[2] = zetanew - zetacur;
    xicur = xinew;
    etacur = etanew;
    zetacur = zetanew;

  } while (!within_tolerance(vector_norm_sq(xidiff, 3), isInElemTolerance_) && ++i < maxNonlinearIter_);

  paraCoords[0] = paraCoords[1] = paraCoords[2] = std::numeric_limits<double>::max();
  double dist = std::numeric_limits<double>::max();

  if (i <maxNonlinearIter_) {
    paraCoords[0] = xinew;
    paraCoords[1] = etanew;
    paraCoords[2] = zetanew;

    std::vector<double> xtmp(3);
    for (int iCom = 0; iCom < 3; iCom++) {
      if (paraCoords[iCom] > 1.0 && paraCoords[iCom] < 1.000001) {
        paraCoords[iCom] = 1.0;
      }
      if (paraCoords[iCom] < -1.0 && paraCoords[iCom] > -1.000001) {
        paraCoords[iCom] = -1.0;
      }
      xtmp[iCom] = paraCoords[iCom];
    }
    dist = parametric_distance(xtmp);
  }
  if (! (dist > 1.0))
    return true;
  else
    return false;
}

/*---------------------
   ~ within_tolerance 
---------------------*/
bool
MapVoxelManager::within_tolerance(const double & val, const double & tol)
{
  return (fabs(val)<tol);
}

/*------------------
  vector_norm_sq 
--------------------*/
double
MapVoxelManager::vector_norm_sq(const double * vect, int len)
{
  double norm_sq = 0.0;
  for (int i = 0; i<len; i++) {
    norm_sq += vect[i] * vect[i];
  }
  return norm_sq;
}

/*----------------------------------
   ~ parametric_distance 
  ----------------------------------*/
double 
MapVoxelManager::parametric_distance(std::vector<double> x)
{
  std::vector<double> y(3);
  for (int i = 0; i<3; ++i) {
    y[i] = fabs(x[i]);
  }

  double d = 0;
  for (int i = 0; i<3; ++i) {
    if (d < y[i]) {
      d = y[i];
    }
  }
  return d;
}

/*--------------------------
   ~ evaluate_temperature
  -------------------------*/
bool
MapVoxelManager::evaluate_temperature(int voxelID, double &interpolatinValue)
{
  double nullMatTemperature = nullMatTemperature_ + 0.00001;
  double lower = caManager_->get_lower_temperature();

  // Double check if whether the cell is with an element
  int tag = elementTagList_[voxelID];
  if (elementTagList_[voxelID] < 0)
  {
    return false;
  }
  int elementID = elementTagList_[voxelID];
  Element * ele = feManager_->elementVector_[elementID];
  if (!ele->birth_) {
    // element is not active  
    return false;
  }
  double * thetaArray = feManager_->thetaArray_;
  bool nonVoid = true;
  for (int i = 0; i < 8; i++)
  {
      double nodalTemp = thetaArray[ele->nID_[i]];
      if (nodalTemp < nullMatTemperature) {
          return false;
      }
  }
  double shapefct[8];
  double thetaIp = 0.0;

  double xicur = parCoordsList_[voxelID * 3];
  double etacur = parCoordsList_[voxelID * 3 + 1];
  double zetacur = parCoordsList_[voxelID * 3 + 2];

  shapefct[0] = 0.125*(1.0 - etacur)*(1.0 - xicur)*(1.0 - zetacur);
  shapefct[1] = 0.125*(1.0 - etacur)*(1.0 + xicur)*(1.0 - zetacur);
  shapefct[2] = 0.125*(1.0 + etacur)*(1.0 + xicur)*(1.0 - zetacur);
  shapefct[3] = 0.125*(1.0 + etacur)*(1.0 - xicur)*(1.0 - zetacur);
  shapefct[4] = 0.125*(1.0 - etacur)*(1.0 - xicur)*(1.0 + zetacur);
  shapefct[5] = 0.125*(1.0 - etacur)*(1.0 + xicur)*(1.0 + zetacur);
  shapefct[6] = 0.125*(1.0 + etacur)*(1.0 + xicur)*(1.0 + zetacur);
  shapefct[7] = 0.125*(1.0 + etacur)*(1.0 - xicur)*(1.0 + zetacur);

    if (ele->birth_) {
      for (int i = 0; i < 8; i++)
      {
        double nodalTemp = thetaArray[ele->nID_[i]];
        thetaIp += shapefct[i] * nodalTemp;
      }

      interpolatinValue = thetaIp;
      return true;
    }
}
