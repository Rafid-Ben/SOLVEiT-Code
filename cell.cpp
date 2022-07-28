#include "cell.h"

/**
 * Creates a cell (node) for use in the tree
 * @param dims - dimensions of the cell
 */
Cell::Cell(vec3<float> dims){
    index = -1;
    numSubcells = 0;
    charge = 0;
    unsignedcharge = 0;
    centerOfCharge = 0;
    dimensions = dims;
    float Ldim[3] = {dimensions.x, dimensions.y, dimensions.z};
    qsort(Ldim, 3, sizeof(float), compareFloat);
    repLength = Ldim[2];
}

/**
 * Deconstructor for cell
 */
Cell::~Cell(){
//    if (firstChild != NULL){
//        delete firstChild;
//        firstChild = NULL;
//    }
//    else if (nextSibling != NULL){
//        delete nextSibling;
//        nextSibling = NULL;
//    }
    if (numSubcells != 0){
        Cell* temp = firstChild;
        Cell* tempNext = temp->nextSibling;
        int i;
        for (i = 0; i < numSubcells - 1; i++){
            delete temp;
            temp = tempNext;
            tempNext = temp->nextSibling;
        }
        delete temp;
    }
}

/**
 * Sets the locations of newly created subcells
 */
void Cell::setLocationOfSubcells(){
    int i;
    vec3<float> subcellDims = dimensions / (double)subcellBlock;
    
    Cell* temp = firstChild;
    int tempIdx;
    
    if (subcellBlock == 2){
        // If 2x2x2 subcells, use faster operations to compute locations
        for (i = 0; i < numSubcells; i++){
            temp->location.x += (i & 1) * subcellDims.x;
            temp->location.y += ((i & 2)>>1) * subcellDims.y;
            temp->location.z += ((i & 4)>>2) * subcellDims.z;
            temp = temp->nextSibling;
        }
    }
    else {
        // If nxnxn subcells, need to determine locations individually
        for (i = 0; i < numSubcells; i++){
            tempIdx = i;
            temp->location.x += (tempIdx % subcellBlock) * subcellDims.x;
            tempIdx /= subcellBlock;
            temp->location.y += (tempIdx % subcellBlock) * subcellDims.y;
            tempIdx /= subcellBlock;
            temp->location.z += (tempIdx % subcellBlock) * subcellDims.z;
            temp = temp->nextSibling;
        }
    }
    
//    firstChild->location = location;
    
//    Cell* temp = firstChild->nextSibling;
//    temp->location.x += subcellDims.x;
//    temp = temp->nextSibling;
//    
//    temp->location.x += subcellDims.x;
//    temp->location.z += subcellDims.z;
//    temp = temp->nextSibling;
//    
//    temp->location.z += subcellDims.z;
//    temp = temp->nextSibling;
//    
//    temp->location.y += subcellDims.y;
//    temp = temp->nextSibling;
//    
//    temp->location.x += subcellDims.x;
//    temp->location.y += subcellDims.y;
//    temp = temp->nextSibling;
//    
//    temp->location += subcellDims;
//    temp = temp->nextSibling;
//    
//    temp->location.y += subcellDims.y;
//    temp->location.z += subcellDims.z;
    
}

/**
 * Creates and initializes subcells when a cell needs splitting
 * @param subcellBlock - number of cells in all dimensions to split a cell into
 * 
 *  Ex: If subcellBlock == 2, cells will be split into 2x2x2 cells
 *      If subcellBlock == 3, cells will be split into 3x3x3 cells
 */
void Cell::generateSubcells(int subcellBlock) {
    this->subcellBlock = subcellBlock;
    
    vec3<float> subcellDims = dimensions / (double)subcellBlock;
    
    numSubcells = subcellBlock * subcellBlock * subcellBlock;
    
    firstChild = new Cell(subcellDims);
    
    Cell* temp = firstChild;
    
    int i;
    for (i = 0; i < (numSubcells - 1); i++){
        temp->location = location;
        temp->nextSibling = new Cell(subcellDims);
        temp = temp->nextSibling;
    }
    temp->location = location;
    temp->nextSibling = nextSibling;
    
    setLocationOfSubcells();
}