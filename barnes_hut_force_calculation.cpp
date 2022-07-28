/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"

/**
 * Construct BarnesHutForceCalculation class and set numThreads and naming
 * @param particleSystem - class containing information about particles
 */
BarnesHutForceCalculation::BarnesHutForceCalculation(ParticleSystem* particleSystem) :
ForceCalculationMethod{particleSystem}
{
    setNumThreads(particleSystem->params.numThreads);
//    switch (particleSystem->params.symType){
//        case (SymmetryType::none):
//            timerName = "Barnes Hut Method";
//            saveName = "_barnesHut.bin";
//            break;
//        case (SymmetryType::rotational):
//            timerName = "Rotational Barnes Hut Method";
//            saveName = "_barnesHutRotational.bin";
//            break;
//        case (SymmetryType::translational):
//            timerName = "Translational Barnes Hut Method";
//            saveName = "_barnesHutTranslational.bin";
//            break;
//    }
}

/**
 * Destructor for BarnesHutForceCalculation
 * Deallocate memory associated with threading
 */
BarnesHutForceCalculation::~BarnesHutForceCalculation() {
//    printf("Derived class destructor called!\n");
    deallocateThreadMemory();
}

/**
 * Allocate memory associated with threading
 */
void BarnesHutForceCalculation::allocateThreadMemory() {
#ifdef _PTHREAD_H
    if (dF != nullptr) delete[] dF;
    dF = new barnesHutForceStruct[numThreads];
#endif
}

/**
 * Deallocate memory associated with threading
 */
void BarnesHutForceCalculation::deallocateThreadMemory() {
#ifdef _PTHREAD_H
    if (dF != nullptr) delete[] dF;
    dF = nullptr;
#endif
}

/**
 * Preset the dimensions of the tree domain used to cover all particles
 * @param boxDimVals - dimensions of the tree domain
 * 
 *  Ex: If boxDimVals = {1, 1, 3}, the tree structure used to cover all
 *      particles will use one cell in the x and y directions and 3 cells in
 *      the z direction. These initial cells will then be subdivided further as
 *      particles are placed in the tree.
 * 
 * Main usage would be testing optimal tree structures
 */
void BarnesHutForceCalculation::setBoxDimVals(vec3<int> boxDimVals){
    forceBoxDims = true;
    forcedBoxDimVals = boxDimVals;
}

/**
 * Set the options for the initial domain type
 * @param toggle - specifies a domain type
 * 
 *  Ex: If toggle == DomainToggle::rectangularPrismCells, the starting domains
 *      will be shaped like rectangular prisms, with length determined by the
 *      size required to cover all particles in the domain.
 *      If toggle == DomainToggle::cubicCells, the starting domains will be
 *      shaped like cubes, and will be stacked in such a way to cover all
 *      particles in the domain.
 * 
 * DomainToggle::cubicCells currently yields better performance
 */
void BarnesHutForceCalculation::setDomainToggle(DomainToggle toggle){
    domainToggle = toggle;
}

/**
 * Set the number of subcells to subdivide a cell into
 * @param subcellBlockSize - number of cells in all dimensions to split a cell into
 * 
 *  Ex: If subcellBlockSize == 2, cells will be split into 2x2x2 cells if necessary
 *      If subcellBlockSize == 3, cells will be split into 3x3x3 cells if necessary
 */
void BarnesHutForceCalculation::setSubcellBlockSize(int subcellBlockSize){
    subcellBlock = subcellBlockSize;
}

/**
 * Set the value of theta (opening angle) for the Barnes-Hut method
 * @param tempTheta - opening angle for the Barnes-Hut method
 */
void BarnesHutForceCalculation::setTheta(float tempTheta){
    theta = tempTheta;
}

/**
 * Toggle whether to ignore outlier particles when creating tree structure
 * @param toggle - true if ignore outlier particles, false if include particles
 * 
 * Benefit likely depends on the particle distribution and location of outliers
 *  If some outlier particles are extremely far away from the majority of 
 *  particles, ignoring the outliers when creating the tree structure has the
 *  potential benefit of a more optimal initial tree structure
 * 
 * Note that all forces will still be calculated, regardless of the option chosen
 */
void BarnesHutForceCalculation::setOutlierToggle(int toggle){
    outlierToggle = toggle;
}

/**
 * Compute effective distance between a particle and a cell (with softening)
 * @param index - index of particle
 * @param cell - pointer to interacting cell
 * @return Effective distance between particle and cell
 */
float BarnesHutForceCalculation::computeEffectiveDistance(int index, Cell* cell){
    float distx = particleSystem->chargedParticles.position[index].x - cell->centerOfCharge.x;
    float disty = particleSystem->chargedParticles.position[index].y - cell->centerOfCharge.y;
    float distz = particleSystem->chargedParticles.position[index].z - cell->centerOfCharge.z;
    
    return sqrt(distx * distx + disty * disty + distz * distz + PhysicalConstant::SOFTENING);
}

/**
 * Locates the subcell to which a particle must be added
 * @param cell - pointer to parent cell
 * @param index - index of particle to be added to subcell
 * @return Integer value indicating the location of subcell
 */
int BarnesHutForceCalculation::locateSubcell(Cell* cell, int index){
    // Determine which subcell to add the body to
    
    if (cell->subcellBlock == 2){
        // With 2x2x2 subcells, can use faster operations to evaluate location
        vec3<float> cellCenter = cell->location + cell->dimensions / 2.0;
        
        int temp = ((particleSystem->chargedParticles.position[index].x > cellCenter.x) + 
                ((particleSystem->chargedParticles.position[index].y > cellCenter.y) << 1) +
                ((particleSystem->chargedParticles.position[index].z > cellCenter.z) << 2));
        
        int temp2;
        if (particleSystem->chargedParticles.position[index].x > cellCenter.x){
            if (particleSystem->chargedParticles.position[index].y > cellCenter.y){
                if (particleSystem->chargedParticles.position[index].z > cellCenter.z){
                    temp2 = 7;
                }
                else
                    temp2 = 3;
            }
            else{
                if (particleSystem->chargedParticles.position[index].z > cellCenter.z){
                    temp2 = 5;
                }
                else
                    temp2 = 1;
            }
        }
        else{
            if (particleSystem->chargedParticles.position[index].y > cellCenter.y){
                if (particleSystem->chargedParticles.position[index].z > cellCenter.z){
                    temp2 = 6;
                }
                else
                    temp2 = 2;
            }
            else{
                if (particleSystem->chargedParticles.position[index].z > cellCenter.z){
                    temp2 = 4;
                }
                else
                    temp2 = 0;
            }
        }
        
        if (temp != temp2){
            printf("error\n");
        }
        
        return ((particleSystem->chargedParticles.position[index].x > cellCenter.x) + 
                ((particleSystem->chargedParticles.position[index].y > cellCenter.y) << 1) +
                ((particleSystem->chargedParticles.position[index].z > cellCenter.z) << 2));
    }
    else {
        // With nxnxn subcells, need to determine location for each dimension
        
        vec3<int> dim_toggles = 0;
        dim_toggles.x = (((double)cell->subcellBlock) / cell->dimensions.x) * 
                (particleSystem->chargedParticles.position[index].x - cell->location.x);
        dim_toggles.y = (((double)cell->subcellBlock) / cell->dimensions.y) * 
                (particleSystem->chargedParticles.position[index].y - cell->location.y);
        dim_toggles.z = (((double)cell->subcellBlock) / cell->dimensions.z) * 
                (particleSystem->chargedParticles.position[index].z - cell->location.z);
        
        return (dim_toggles.x + cell->subcellBlock * dim_toggles.y + 
                cell->subcellBlock * cell->subcellBlock * dim_toggles.z);
    }
}

/**
 * Adds a particle to a cell. If a particle already exists in the cell, the cell
 * is sub-divided and both particles are added to subcells
 * @param cell - pointer to cell
 * @param index - index of particle to add to cell
 */
void BarnesHutForceCalculation::addToCell(Cell* cell, int index){
    if (cell->index == -1){
        cell->index = index;
    }
    else{
        cell->generateSubcells(subcellBlock);
        
        // The current cell's body must now be re-added to one of its subcells
        int subcellIndex1 = locateSubcell(cell, cell->index);
        Cell* temp = cell->firstChild;
        for (int i = 0; i < subcellIndex1; i++){
            temp = temp->nextSibling;
        }
        temp->index = cell->index;
        
        // Locate subcell for new body
        int subcellIndex2 = locateSubcell(cell, index);
        
        if (subcellIndex1 == subcellIndex2){
            addToCell(temp, index);
        }
        else{
            temp = cell->firstChild;
            for (int i = 0; i < subcellIndex2; i++){
                temp = temp->nextSibling;
            }
            temp->index = index;
        }
    }
}

/**
 * Generates the tree for the entire system of particles
 */
void BarnesHutForceCalculation::generateTree(){
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);
    
    int i, ii, j, k;
    
    // Set locations of the largest cells
    rootCell->location = roots;
    Cell* temp = rootCell;
    for (i = 0; i < domainDimVals.x; i++){
        for (j = 0; j < domainDimVals.y; j++){
            for (k = 0; k < domainDimVals.z; k++){
                if (i != 0 || j != 0 || k != 0){
                    temp = temp->nextSibling;
                    temp->location.x = roots.x + i * dims.x;
                    temp->location.y = roots.y + j * dims.y;
                    temp->location.z = roots.z + k * dims.z;
                }
            }
        }
    }
    
    // Add all particles to a cell in the tree structure
    vec3<int> domainValTemp;
    int domainFound;
    Cell* cell;
    for (ii = firstParticle; ii < lastParticle; ii++){
        if (particleSystem->chargedParticles.acceleration[ii].w == inboundCheck ||
                particleSystem->chargedParticles.acceleration[ii].w == affectsAccelerationCheck){
            
            // Find highest level cell
            domainValTemp.x = (particleSystem->chargedParticles.position[ii].x - roots.x) / dims.x;
            domainValTemp.y = (particleSystem->chargedParticles.position[ii].y - roots.y) / dims.y;
            domainValTemp.z = (particleSystem->chargedParticles.position[ii].z - roots.z) / dims.z;
            cell = rootCell;
            domainFound = 0;

            i = 0;
            while (i < domainDimVals.x && domainFound == 0){
                j = 0;
                while (j < domainDimVals.y && domainFound == 0){
                    k = 0;
                    while (k < domainDimVals.z && domainFound == 0){
                        if (i != 0 || j != 0 || k != 0){
                            cell = cell->nextSibling;
                        }
                        if (i == domainValTemp.x && j == domainValTemp.y && k == domainValTemp.z){
                            domainFound = 1;
                        }
                        k++;
                    }
                    j++;
                }
                i++;
            }

            if (domainFound == 0){
                printf("Error: domain not found!\n");
            }

            // Add particle to cell
            while (cell->numSubcells != 0) {
                int subcellIndex = locateSubcell(cell, ii);
                temp = cell->firstChild;
                for (j = 0; j < subcellIndex; j++){
                    temp = temp->nextSibling;
                }
                cell = temp;
            }

            addToCell(cell, ii);
        }
    }
}

/**
 * Computes the total charge and the center of charge of the current cell
 * @param cell - pointer to cell
 */
void BarnesHutForceCalculation::computeCellProperties(Cell* cell){
    if (cell->index == 0){
        //printf("cell.index = %d\n", cell.index);
    }
    
    if (cell->numSubcells == 0){
        if (cell->index != -1){
            // Just a single particle in a single cell
//            cell->charge = getCharge(particles.velocity[cell->index].w, 0.0f);
//            cell->unsignedcharge = abs(getCharge(particles.velocity[cell->index].w, 0.0f));
            cell->charge = particleSystem->chargedParticles.getCharge(cell->index);
            cell->unsignedcharge = std::abs(particleSystem->chargedParticles.getCharge(cell->index));
            cell->centerOfCharge.x = particleSystem->chargedParticles.position[cell->index].x;
            cell->centerOfCharge.y = particleSystem->chargedParticles.position[cell->index].y;
            cell->centerOfCharge.z = particleSystem->chargedParticles.position[cell->index].z;
        }
    }
    else{
        int i;
        vec3<float> t(0,0,0);
        Cell* temp = cell->firstChild;
        for (i = 0; i < cell->numSubcells; i++){
            // Compute cell properties for each subcell
            computeCellProperties(temp);
//            if (cell.subcells[i] != NULL) {
            cell->charge += temp->charge;
            cell->unsignedcharge += temp->unsignedcharge;

            t += temp->centerOfCharge * std::abs(temp->charge);
//            }
            
//            if (cell->unsignedcharge != 0){
//                cell->centerOfCharge = t / cell->unsignedcharge;
//            }
//            else{
//                cell->centerOfCharge = cell->location + cell->dimensions / 2;
//            }
            
            temp = temp->nextSibling;
        }
        
        // Set the location of the center of charge
        if (cell->unsignedcharge != 0){
            cell->centerOfCharge = t / cell->unsignedcharge;
        }
        else{
            cell->centerOfCharge = cell->location + cell->dimensions / 2.0;
        }
    }
    
    vec3<float> center = cell->centerOfCharge - (cell->location + cell->dimensions / 2.0);
        
    cell->delta = sqrt(center.x * center.x + center.y * center.y + center.z * center.z);
}

/**
 * Computes the force between a particle and a cell
 * @param cell - pointer to cell
 * @param index - index of particle
 * @param accel - reference to acceleration of particle
 * @param deff - distance between particle and cell
 */
inline void BarnesHutForceCalculation::computeForceFromCell(Cell* cell, int index, vec3<float>& accel, float deff){
//    float deff = computeEffectiveDistance(index, cell);
    
    // Compute electric force according to Coulomb's law
    float f = cell->charge / (deff * deff * deff);
    
    accel.x += f * (particleSystem->chargedParticles.position[index].x - cell->centerOfCharge.x);
    accel.y += f * (particleSystem->chargedParticles.position[index].y - cell->centerOfCharge.y);
    accel.z += f * (particleSystem->chargedParticles.position[index].z - cell->centerOfCharge.z);
}

/**
 * Computes the force between the particles in the system, using the clustering
 *  approximation for long distant forces
 * @param cell - pointer to cell
 * @param index - index of particle
 * @param accel - reference to acceleration of particle
 */
void BarnesHutForceCalculation::computeForceFromTree(Cell* cell, int index, vec3<float>& accel){
    
    Cell* temp = cell;
    float d;
    while (temp != NULL){
        if (temp->firstChild == NULL){
            if (temp->index != -1 && temp->index != index){
                d = computeEffectiveDistance(index, temp);
//                if (index == 28){
//                    printf("Only one particle: %p, dist: %.9g, charge: %.9g\n",(void*)temp, d,temp->charge);
//                }
                computeForceFromCell(temp, index, accel, d);
            }
            temp = temp->nextSibling;
        }
        else {
            d = computeEffectiveDistance(index, temp);
//            if (d > ((temp->repLength / theta) + temp->delta)){
            if ((theta * (d - temp->delta)) > temp->repLength){
//                if (index == 28){
//                    printf("Far enough away: %p, dist: %.9g, charge: %.9g\n",(void*)temp, d,temp->charge);
//                }
                computeForceFromCell(temp, index, accel, d);
                temp = temp->nextSibling;
            }
            else {
//                if (index == 28){
//                    printf("Too close: %p, dist: %.9g, charge: %.9g\n",(void*)temp, d,temp->charge);
//                }
                temp = temp->firstChild;
            }
        }
    }
}

/**
 * Computes the boundaries of the simulation, and sets the dimensions for the
 *  highest level of the cell structure
 */
void BarnesHutForceCalculation::computeWorldLimits(){
    
    float rootBoxSize;
    vec3<float> minOutliers, maxOutliers;
    vec3<int> boxDimVals, dimInd;
    int avgPPB;
    DomainBoundsFunction::computeDomainBounds(particleSystem->chargedParticles,
            outlierCheck, boxDimVals, roots, minOutliers, maxOutliers, dims,
            numOutliers, forceBoxDims, forcedBoxDimVals, domainToggle, avgPPB, dimInd);

    float subtr = sqrt(PhysicalConstant::SOFTENING);
    
//    while (roots.x > minOutliers.x && std::abs(roots.x - minOutliers.x) > subtr){
    while (roots.x > minOutliers.x){
        roots.x -= dims.x;
        boxDimVals[0]++;
    }
//    while (roots.y > minOutliers.y && std::abs(roots.y - minOutliers.y) > subtr){
    while (roots.y > minOutliers.y){
        roots.y -= dims.y;
        boxDimVals[1]++;
    }
//    while (roots.z > minOutliers.z && std::abs(roots.z - minOutliers.z) > subtr){
    while (roots.z > minOutliers.z){
        roots.z -= dims.z;
        boxDimVals[2]++;
    }
//    while (roots.x + boxDimVals[0] * dims.x < maxOutliers.x && std::abs(roots.x + boxDimVals[0] * dims.x - maxOutliers.x) > subtr){
    while (roots.x + boxDimVals[0] * dims.x < maxOutliers.x){
        boxDimVals[0]++;
    }
//    while (roots.y + boxDimVals[1] * dims.y < maxOutliers.y && std::abs(roots.y + boxDimVals[1] * dims.y - maxOutliers.y) > subtr){
    while (roots.y + boxDimVals[1] * dims.y < maxOutliers.y){
        boxDimVals[1]++;
    }
//    while (roots.z + boxDimVals[2] * dims.z < maxOutliers.z && std::abs(roots.z + boxDimVals[2] * dims.z - maxOutliers.z) > subtr){
    while (roots.z + boxDimVals[2] * dims.z < maxOutliers.z){
        boxDimVals[2]++;
    }

    if (domainToggle == DomainToggle::cubicCells){
        rootBoxSize = dims.x; // Should equal dims.y and dims.z as well
        int domainCheck = std::min(boxDimVals[0], std::min(boxDimVals[1], boxDimVals[2]));
    //    printf("dim.x = %d, .y = %d, .z = %d\n",boxDimVals[0], boxDimVals[1], boxDimVals[2]);
        while (domainCheck / 2 > 0){
            domainCheck = (domainCheck + 1) / 2;
            rootBoxSize *= 2;
            boxDimVals[0] = (boxDimVals[0] + 1) / 2;
            boxDimVals[1] = (boxDimVals[1] + 1) / 2;
            boxDimVals[2] = (boxDimVals[2] + 1) / 2;
        }
        dims = rootBoxSize;
    }

    domainDimVals.x = boxDimVals[0];
    domainDimVals.y = boxDimVals[1];
    domainDimVals.z = boxDimVals[2];
    
//    if (domainDimVals.x * domainDimVals.y * domainDimVals.z > 1){
//        printf("test\n");
//    }
//    printf("dim.x = %d, .y = %d, .z = %d\n",domainDimVals.x, domainDimVals.y, domainDimVals.z);
//    printf("dims.x = %.9g, .y = %.9g, .z = %.9g\n",dims.x,dims.y,dims.z);
//    if (domainDimVals.z == 5){
////        printf("test\n");
//    }
}

/**
 * Multi-threaded function to call the Barnes-Hut force computation method
 * @param dF - struct that holds a pointer to the class and the thread number
 * @return nullptr
 */
void* BarnesHutForceCalculation::barnesHutKernel(void* dF){
    barnesHutForceStruct dfc = *((barnesHutForceStruct*)dF);
    dfc.d->barnesHutKernelCalculation(dfc.threadNum);
    return nullptr;
}

/**
 * Computes the force on all particles using the Barnes-Hut approximation
 * @param threadNum - thread number for multi-threaded versions
 */
void BarnesHutForceCalculation::barnesHutKernelCalculation(int threadNum){
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    vec3<float> accel;
    double accel_const;
    int i;
    for (i = rangevalues[threadNum].x; i < rangevalues[threadNum].y; i++){
//        if (rangevalues[threadNum].y == 485 && i == 159){
//            printf("test\n");
//        }
//        if (i == 28){
//            printf("test\n");
//        }
        if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
            accel = 0.0f;
            computeForceFromTree(rootCell, i, accel);

            accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / 
                    particleSystem->chargedParticles.getMass(i);
            accel *= accel_const;

            particleSystem->chargedParticles.acceleration[i].x += accel.x;
            particleSystem->chargedParticles.acceleration[i].y += accel.y;
            particleSystem->chargedParticles.acceleration[i].z += accel.z;
        }
        else {
            particleSystem->chargedParticles.acceleration[i].x += 0;
            particleSystem->chargedParticles.acceleration[i].y += 0;
            particleSystem->chargedParticles.acceleration[i].z += 0;
        }
    }
}

/**
 * Updates the acceleration of all particles due to particle-particle interactions
 */
void BarnesHutForceCalculation::updateAcceleration(){
    setTimerName();
    
    if (particleSystem == nullptr){
        printf("Particle system not initialized.\n");
        exit(1);
    }
    
    int numParticles = lastParticle - firstParticle;
    int i, j, k;
    
    Timer precalcTimer, accelTimer;
    float precalcTime, accelTime;
    
    precalcTimer.startTimer();
    
    outlierCheck = new int[numParticles];
    
    if (outlierToggle == 1){
        DomainBoundsFunction::isOutlier(particleSystem->chargedParticles, 
                outlierCheck, numOutliers);
    }
    else if (outlierToggle == 0){
        numOutliers = 0;
        for (i = 0; i < numParticles; i++){
            outlierCheck[i] = 0;
        }
    }
    
    computeWorldLimits();
    
    delete[] outlierCheck;
    
    rootCell = new Cell(dims);
    Cell* temp = rootCell;
    
    for (i = 0; i < domainDimVals.x; i++){
        for (j = 0; j < domainDimVals.y; j++){
            for (k = 0; k < domainDimVals.z; k++){
                if (i != 0 || j != 0 || k != 0){
                    temp->nextSibling = new Cell(dims);
                    temp = temp->nextSibling;
                }
            }
        }
    }
    
    generateTree();
    
    computeCellProperties(rootCell);
    temp = rootCell;
    for (i = 0; i < domainDimVals.x; i++){
        for (j = 0; j < domainDimVals.y; j++){
            for (k = 0; k < domainDimVals.z; k++){
                if (i != 0 || j != 0 || k != 0){
                    temp = temp->nextSibling;
                    computeCellProperties(temp);
                }
            }
        }
    }
    
    precalcTimer.pauseTimer();
    accelTimer.startTimer();
    
#ifdef _PTHREAD_H
//    int numParticles = lastParticle - firstParticle;
    for (k = 0; k < numThreads; k++){
        dF[k].d = this;
        dF[k].threadNum = k;
        rangevalues[k].x = (k * numParticles) / numThreads + firstParticle;
        rangevalues[k].y = ((k + 1) * numParticles) / numThreads + firstParticle;
        rangevalues[k].z = firstParticle;
        rangevalues[k].w = lastParticle;
        
        switch (particleSystem->params.symType){
            case (SymmetryType::none):
                pthread_create(&threads[k], nullptr, barnesHutKernel, &dF[k]);
                break;
            case (SymmetryType::rotational):
                pthread_create(&threads[k], NULL, rotationalBarnesHutKernel, &dF[k]);
                break;
            case (SymmetryType::translational):
                pthread_create(&threads[k], NULL, translationalBarnesHutKernel, &dF[k]);
                break;
        }
    }

    for (k = 0; k < numThreads; k++) {
        pthread_join(threads[k], NULL);
    }
#else
    rangevalues[0].x = firstParticle;
    rangevalues[0].y = lastParticle;
    rangevalues[0].z = firstParticle;
    rangevalues[0].w = lastParticle;
    switch (particleSystem->params.symType){
        case (SymmetryType::none):
            barnesHutKernelCalculation();
            break;
        case (SymmetryType::rotational):
            rotationalBarnesHutKernelCalculation();
            break;
        case (SymmetryType::translational):
            translationalBarnesHutKernelCalculation();
            break;
    }
#endif
    
    accelTime = accelTimer.getTimer();
    precalcTimer.startTimer();
    
    /*Note use delete, not delete[] since rootCell is not an array of objects*/
    temp = rootCell->nextSibling;
    Cell* tempNext;
    if (temp != NULL){
        tempNext = temp->nextSibling;
    }
    delete rootCell;
    for (i = 0; i < domainDimVals.x; i++){
        for (j = 0; j < domainDimVals.y; j++){
            for (k = 0; k < domainDimVals.z; k++){
                if (i != 0 || j != 0 || k != 0){
                    delete temp;
                    temp = tempNext;
                    if (temp != NULL){
                        tempNext = temp->nextSibling;
                    }
                }
            }
        }
    }
    
    precalcTime = precalcTimer.getTimer();
    
    if (timeInfoFile != nullptr){
        fwrite(&numParticles, sizeof(*&numParticles), 1, timeInfoFile);
        fwrite(&precalcTime, sizeof(*&precalcTime), 1, timeInfoFile);
        fwrite(&accelTime, sizeof(*&accelTime), 1, timeInfoFile);
    }
    
}
