/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"

void MultipoleMethodForceCalculation::updateAccelerationTranslational(){
    int numParticles = lastParticle - firstParticle;
    int i, numLevel, numBoxIndex, numBoxIndexOld;
    
    if (timeInfoFile != nullptr){
        fwrite(&numParticles, sizeof(*&numParticles), 1, timeInfoFile);
    }
    
    timing.resetTiming();
    
    timing.total_timer.startTimer();
    timing.timer.startTimer();
    
    numBoxIndexLeaf = 0;
    numBoxIndexTotal = 0;
    numBoxIndexFull = 0;
    maxLevel = 0;
    numLevelSC = 0;
    numDomains = 0;
    numOutlierDomains = 0;
    numOutliers = 0;
    outlierNumBoxIndexLeaf = 0;
    
    variables->allocateMMvariables(particleSystem->chargedParticles.maxNumParticles, 
            numBoxIndexLeaf, numBoxIndexTotal, numBoxIndexFull, numExpansions,
            maxLevel, numLevelSC, numDomains, numOutlierDomains);
    
    if (outlierToggle == 1){
        DomainBoundsFunction::isOutlier(particleSystem->chargedParticles,
                variables->outlierCheck, numOutliers);
    }
    else if (outlierToggle == 0){
        numOutliers = 0;
        for (i = 0; i < numParticles; i++){
            variables->outlierCheck[i] = 0;
        }
    }

    variables->allocateOutlierVariables(numOutlierDomains, numOutliers, outlierNumBoxIndexLeaf);
    
    setDomainSizeTranslations();
    
    std::sort(variables->outlierList, variables->outlierList + numOutliers, compareOutlierInfo);
    
    numOutlierDomains = 0;
    for (i = 0; i < numOutliers; i++) {
        if (i == 0 || variables->outlierList[i].signedMorton != variables->outlierList[i - 1].signedMorton) {
            numOutlierDomains++;
        }
    }
    variables->allocateOutlierVariables(numOutlierDomains, numOutliers, outlierNumBoxIndexLeaf);

    setOptimumLevelTranslations();
    
    if (timeInfoFile != nullptr){
        fwrite(&numOutliers, sizeof(*&numOutliers), 1, timeInfoFile);
        fwrite(&boxDimVals[0], sizeof(*&boxDimVals[0]), 1, timeInfoFile);
        fwrite(&boxDimVals[1], sizeof(*&boxDimVals[1]), 1, timeInfoFile);
        fwrite(&boxDimVals[2], sizeof(*&boxDimVals[2]), 1, timeInfoFile);
        fwrite(&maxLevel, sizeof(*&maxLevel), 1, timeInfoFile);
    }

    if (boxInfo != nullptr){
        fwrite(&boxMin.x, sizeof (*&boxMin.x), 1, boxInfo);
        fwrite(&boxMin.y, sizeof (*&boxMin.y), 1, boxInfo);
        fwrite(&boxMin.z, sizeof (*&boxMin.z), 1, boxInfo);
        fwrite(&rootBoxSize, sizeof (*&rootBoxSize), 1, boxInfo);
        fwrite(&boxDimVals[0], sizeof (*&boxDimVals[0]), 1, boxInfo);
        fwrite(&boxDimVals[1], sizeof (*&boxDimVals[1]), 1, boxInfo);
        fwrite(&boxDimVals[2], sizeof (*&boxDimVals[2]), 1, boxInfo);
        fwrite(&maxLevel, sizeof (*&maxLevel), 1, boxInfo);
    }
    
    if (maxLevel == -1) {
        
        timing.setup_time = timing.timer.getTimer();
        timing.timer.resetTimer();
        
        int k;
        if (vectorizedToggle){
            int tempNumParticles;
            DirectForceFunction::setJpdata(tempNumParticles, particleSystem->chargedParticles, jptc);
            for (k = 0; k < numThreads; k++){
                rangevalues[k].z = 0;
                rangevalues[k].w = tempNumParticles;
            }
        }
        else {
            for (k = 0; k < numThreads; k++){
                rangevalues[k].z = firstParticle;
                rangevalues[k].w = lastParticle;
            }
        }
        
        // Small number of particles -- use direct method
        direct();
    } else {
        variables->allocateMMvariables(particleSystem->chargedParticles.maxNumParticles, 
                numBoxIndexLeaf, numBoxIndexTotal, numBoxIndexFull, numExpansions,
                maxLevel, numLevelSC, numDomains, numOutlierDomains);
        
        sortParticlesTranslations();
        
        if (vectorizedToggle){
            int tempNumParticles;
            DirectForceFunction::setJpdata(tempNumParticles, particleSystem->chargedParticles, jptc);
        }
        
        countNonEmptyBoxesTranslations();
        
        variables->allocateMMvariables(particleSystem->chargedParticles.maxNumParticles, 
                numBoxIndexLeaf, numBoxIndexTotal, numBoxIndexFull, numExpansions,
                maxLevel, numLevelSC, numDomains, numOutlierDomains);

        numLevel = maxLevel;

        for (i = firstParticle; i < lastParticle; i++) {
            particleSystem->chargedParticles.acceleration[i].x = 0;
            particleSystem->chargedParticles.acceleration[i].y = 0;
            particleSystem->chargedParticles.acceleration[i].z = 0;
        }
        
        timing.setup_time = timing.timer.getTimer();
        timing.timer.resetTimer();
        
        if (numLevel > 0){
            variables->levelOffset[numLevel - 1] = 0;
            
            getBoxDataTranslations(numBoxIndex, 0);
            
            // P2P
            getInteractionListTranslations(numBoxIndex, numLevel, FmmInteractionType::P2P);

            if (numOutliers > 0) {
                outlier_p2p();
            }
            
            p2p(numBoxIndex, 0);
            
            // Interaction list for interactions between outliers and particles close by
            if (numOutliers > 0) {
                getInteractionListTranslations(outlierNumBoxIndex, numLevel, FmmInteractionType::P2P_Outliers);
                p2outlier_p(outlierNumBoxIndex, 0);
            }
            
            // P2M
            p2m(numBoxIndex, 0);
        } else {
            variables->levelOffset[0] = 0;
            getBoxDataTranslations(numBoxIndex, 1);
            // P2P for supercells
            getInteractionListTranslations(numBoxIndex, numLevel, FmmInteractionType::P2P_Supercell);
            if (numOutliers > 0) {
                outlier_p2p();
            }

            p2p(numBoxIndex, 1);

            // Interaction list for interactions between outliers and particles close by
            if (numOutliers > 0) {
                getInteractionListTranslations(outlierNumBoxIndex, numLevel, FmmInteractionType::P2P_SupercellOutliers);
                p2outlier_p(outlierNumBoxIndex, 1);
            }

            // P2M for supercells

            p2m(numBoxIndex, 1);
        }
        
        numLevel = maxLevel;

        int tempLevel = maxLevel;
        int tempToggle = (maxLevel == 0) ? 1 : 0;
        
        if (maxLevel + numLevelSC <= 0) {
            getBoxIndexMaskTranslations(numBoxIndex, numLevel, 0);
        }
        
        for (numLevel = maxLevel - 1; numLevel >= std::max(0, 2 - numLevelSC); numLevel--){
            if (treeOrFMM == 0) {

                // M2P at lower levels

                getInteractionListTranslations(numBoxIndex, numLevel + 1, FmmInteractionType::M2L_Subcell);

                m2p(numBoxIndex, numLevel + 1, 0);
            }
            
            if (numOutliers > 0) {
                getInteractionListTranslations(outlierNumBoxIndex, numLevel + 1, FmmInteractionType::M2L_SubcellOutliers);
                m2outlier_p(outlierNumBoxIndex, numLevel + 1, 0);
            }

            numBoxIndexOld = numBoxIndex;

            if (numLevel == 0) {
                getBoxDataOfParentTranslations(numBoxIndex, numLevel, 1);
            } else {
                getBoxDataOfParentTranslations(numBoxIndex, numLevel, 0);
            }

            m2m(numBoxIndex, numBoxIndexOld, numLevel, 0);

            tempLevel = numLevel;
            tempToggle = 0;
        }
        
        for (numLevel = 1; numLevel <= numLevelSC - 1; numLevel++) {
            if (treeOrFMM == 0) {
                // M2P at supercell levels

                getInteractionListTranslations(numBoxIndex, numLevel - 1, FmmInteractionType::M2L_Supercell);

                m2p(numBoxIndex, numLevel - 1, 1);
            }
            if (numOutliers > 0) {
                getInteractionListTranslations(outlierNumBoxIndex, numLevel - 1, FmmInteractionType::M2L_SupercellOutliers);
                m2outlier_p(outlierNumBoxIndex, numLevel - 1, 1);
            }

            numBoxIndexOld = numBoxIndex;

            getBoxDataOfParentTranslations(numBoxIndex, numLevel, 1);

            m2m(numBoxIndex, numBoxIndexOld, numLevel, 1);

            tempLevel = numLevel;
            tempToggle = 1;
        }
        
        numLevel = tempLevel;
        FmmInteractionType tempInteractionType;

        if (numOutliers > 0) {
            if (numLevel == 0) tempToggle = 1;
            if (tempToggle == 0){
//                tempInteractionType = FmmInteractionType::M2L_SubcellTopLevelOutliers;
                tempInteractionType = FmmInteractionType::M2L_SubcellOutliers;
            }
            else{
//                tempInteractionType = FmmInteractionType::M2L_SupercellTopLevelOutliers;
                tempInteractionType = FmmInteractionType::M2L_SupercellOutliers;
            }
            getInteractionListTranslations(outlierNumBoxIndex, tempLevel, tempInteractionType); // Using outlierNumBoxIndex instead of numBoxIndex is correct
            m2outlier_p(outlierNumBoxIndex, tempLevel, tempToggle);
        }
        
        if (treeOrFMM == 0) {
            if (numLevel == 0) tempToggle = 1;
            if (tempToggle == 0){
//                tempInteractionType = FmmInteractionType::M2L_SubcellTopLevel;
                tempInteractionType = FmmInteractionType::M2L_Subcell;
            }
            else{
//                tempInteractionType = FmmInteractionType::M2L_SupercellTopLevel;
                tempInteractionType = FmmInteractionType::M2L_Supercell;
            }
            getInteractionListTranslations(numBoxIndex, tempLevel, tempInteractionType);

            m2p(numBoxIndex, tempLevel, tempToggle);
            
        } else if (treeOrFMM == 1) {
            if (numLevel == 0) tempToggle = 1;
            if (tempToggle == 0){
//                tempInteractionType = FmmInteractionType::M2L_SubcellTopLevel;
                tempInteractionType = FmmInteractionType::M2L_Subcell;
            }
            else{
//                tempInteractionType = FmmInteractionType::M2L_SupercellTopLevel;
                tempInteractionType = FmmInteractionType::M2L_Supercell;
            }
            getInteractionListTranslations(numBoxIndex, numLevel, tempInteractionType);

            m2l(numBoxIndex, numLevel, tempToggle, 1);

            // L2L

            if (maxLevel + numLevelSC >= 2) {

                for (numLevel = numLevelSC - 2; numLevel >= 0; numLevel--) {
                    numBoxIndex = variables->levelOffset[maxLevel + numLevel + 1] - variables->levelOffset[maxLevel + numLevel];

                    l2l(numBoxIndex, numLevel, 1);

                    getBoxIndexMaskTranslations(numBoxIndex, numLevel, 1);

                    getInteractionListTranslations(numBoxIndex, numLevel, FmmInteractionType::M2L_Supercell);

                    m2l(numBoxIndex, numLevel, 1);
                }

                for (numLevel = std::max(1, 3 - numLevelSC); numLevel <= maxLevel; numLevel++) {

                    if (numLevel == 1) numBoxIndex = variables->levelOffset[maxLevel] - variables->levelOffset[0];
                    else numBoxIndex = variables->levelOffset[numLevel - 2] - variables->levelOffset[numLevel - 1];

                    l2l(numBoxIndex, numLevel, 0);

                    getBoxIndexMaskTranslations(numBoxIndex, numLevel, 0);

                    // M2L at lower levels

                    getInteractionListTranslations(numBoxIndex, numLevel, FmmInteractionType::M2L_Subcell);

                    m2l(numBoxIndex, numLevel, 0);
                }

            }

            // L2P
            if (maxLevel == 0) {
                l2p(numBoxIndex, 1);
            } else {
                l2p(numBoxIndex, 0);
            }
        }

        unsortParticles();
    }
    
    timing.total_time = timing.total_timer.getTimer();
    timing.writeTiming(timeInfoFile);
}

void MultipoleMethodForceCalculation::setDomainSizeTranslations(){
    int i, j, nx, ny, nz;
    vec3<float> minOutliers, maxOutliers, dims;
    vec3<int> dimInd;
    int avgPPB;
    DomainToggle domainToggle = DomainToggle::cubicCells;
    vec3<int> tempBoxDimVals;
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);
    
    dims = particleSystem->params.domainSize;
    
    DomainBoundsFunction::computeDomainBoundsTranslations(particleSystem->chargedParticles, 
            variables->outlierCheck, tempBoxDimVals, boxMin, minOutliers, maxOutliers,
            dims, numOutliers, forceBoxDims, forcedBoxDimVals, domainToggle,
            avgPPB, dimInd, particleSystem->params.symmetryAxes);
    
    rootBoxSize = dims.x; // dims should all be equal
    
    boxDimVals[0] = tempBoxDimVals[0];
    boxDimVals[1] = tempBoxDimVals[1];
    boxDimVals[2] = tempBoxDimVals[2];
    numDomains = boxDimVals[0] * boxDimVals[1] * boxDimVals[2];
    
    // Not sure if this would ever lead to rounding errors, so below is safer
    //numLevelSC = ceil(log(boxDimVals[dimInd[2]]) / log(3)) + 1;

    int powThree = 1;
    numLevelSC = 1;
    while (boxDimVals[dimInd[2]] > powThree) {
        numLevelSC++;
        powThree *= 3;
    }

    variables->allocateDomainVariables(numLevelSC, numDomains, avgPPB);
    variables->numDomainsLevel[0] = numDomains;
    float boxSizeSC = rootBoxSize;
    int temp = 3;
    variables->dimValsSC[0].x = boxDimVals[0];
    variables->dimValsSC[0].y = boxDimVals[1];
    variables->dimValsSC[0].z = boxDimVals[2];
    variables->boxMinSC[0].x = boxMin.x;
    variables->boxMinSC[0].y = boxMin.y;
    variables->boxMinSC[0].z = boxMin.z;

    variables->boxMinIndexSC[0].x = 0;
    variables->boxMinIndexSC[0].y = 0;
    variables->boxMinIndexSC[0].z = 0;
    variables->offsetSC[0].x = (variables->dimValsSC[0].x % 3) % 2;
    variables->offsetSC[0].y = (variables->dimValsSC[0].y % 3) % 2;
    variables->offsetSC[0].z = (variables->dimValsSC[0].z % 3) % 2;

    // Calculate domain offsets and dimensions for higher level supercells
    for (i = 1; i < numLevelSC; i++) {
        variables->dimValsSC[i].x = (variables->dimValsSC[i - 1].x + 2) / 3;
        variables->dimValsSC[i].y = (variables->dimValsSC[i - 1].y + 2) / 3;
        variables->dimValsSC[i].z = (variables->dimValsSC[i - 1].z + 2) / 3;
        variables->numDomainsLevel[i] = variables->dimValsSC[i].x * variables->dimValsSC[i].y * variables->dimValsSC[i].z;
        variables->offsetSC[i].x = (variables->dimValsSC[i].x % 3) % 2;
        variables->offsetSC[i].y = (variables->dimValsSC[i].y % 3) % 2;
        variables->offsetSC[i].z = (variables->dimValsSC[i].z % 3) % 2;
        variables->boxMinSC[i].x = variables->boxMinSC[i - 1].x - variables->offsetSC[i - 1].x * boxSizeSC;
        variables->boxMinSC[i].y = variables->boxMinSC[i - 1].y - variables->offsetSC[i - 1].y * boxSizeSC;
        variables->boxMinSC[i].z = variables->boxMinSC[i - 1].z - variables->offsetSC[i - 1].z * boxSizeSC;
        boxSizeSC *= 3;

        variables->boxMinIndexSC[i].x = variables->boxMinIndexSC[i - 1].x - temp * variables->offsetSC[i].x;
        variables->boxMinIndexSC[i].y = variables->boxMinIndexSC[i - 1].y - temp * variables->offsetSC[i].y;
        variables->boxMinIndexSC[i].z = variables->boxMinIndexSC[i - 1].z - temp * variables->offsetSC[i].z;
        temp *= 3;
    }
    
    createDomainMortonLink();
    
    int tempIndex;

    vec3<int> tempDomainIndex = {-1, -1, -1};
    int tempMorton3;

    numOutliers = 0;
    unsigned int tempSignedMorton3;
    vec3<int> tempOutlierDomainIndex = variables->offsetSC[0];
    float tempPos;
    
    int numParticles = lastParticle - firstParticle;

    // Assign each particle to a certain domain
    for (i = 0; i < numParticles; i++) {
//        if (i == 113){
//            printf("test\n");
//        }
        if (particleSystem->chargedParticles.acceleration[i + firstParticle].w == inboundCheck ||
                particleSystem->chargedParticles.acceleration[i + firstParticle].w == affectsAccelerationCheck){
        
            tempPos = (particleSystem->chargedParticles.position[i + firstParticle].x - boxMin.x) / rootBoxSize;
            nx = int(tempPos);
            if (tempPos < 0) {
                nx--;
            }
            tempPos = (particleSystem->chargedParticles.position[i + firstParticle].y - boxMin.y) / rootBoxSize;
            ny = int(tempPos);
            if (tempPos < 0) {
                ny--;
            }
            tempPos = (particleSystem->chargedParticles.position[i + firstParticle].z - boxMin.z) / rootBoxSize;
            nz = int(tempPos);
            if (tempPos < 0) {
                nz--;
            }
            if (rootBoxSize == 0){
                nx = 0;
                ny = 0;
                nz = 0;
            }

            //        nx = int((particles.pos[i + initial].x - boxMin.x) / rootBoxSize);
            //        ny = int((particles.pos[i + initial].y - boxMin.y) / rootBoxSize);
            //        nz = int((particles.pos[i + initial].z - boxMin.z) / rootBoxSize);

            if (nx >= boxDimVals[0] || nx < 0 || ny >= boxDimVals[1] || ny < 0 || nz >= boxDimVals[2] || nz < 0) {
                nx += variables->offsetSC[0].x;
                ny += variables->offsetSC[0].y;
                nz += variables->offsetSC[0].z;
                if (nx != tempOutlierDomainIndex.x || ny != tempOutlierDomainIndex.y || nz != tempOutlierDomainIndex.z) {
                    tempOutlierDomainIndex.x = nx;
                    tempOutlierDomainIndex.y = ny;
                    tempOutlierDomainIndex.z = nz;
                    MultipoleMethodFunction::signedMorton3(tempOutlierDomainIndex, tempSignedMorton3);
                }
                variables->outlierList[numOutliers].index = i + firstParticle;
                variables->outlierList[numOutliers].signedMorton = tempSignedMorton3;
                numOutliers++;
            } else {
                if (nx >= boxDimVals[0]) nx--;
                if (ny >= boxDimVals[1]) ny--;
                if (nz >= boxDimVals[2]) nz--;
                nx += variables->offsetSC[0].x;
                ny += variables->offsetSC[0].y;
                nz += variables->offsetSC[0].z;
                if (nx != tempDomainIndex.x || ny != tempDomainIndex.y || nz != tempDomainIndex.z) {
                    tempDomainIndex.x = nx;
                    tempDomainIndex.y = ny;
                    tempDomainIndex.z = nz;
                    MultipoleMethodFunction::morton3(tempDomainIndex, tempMorton3, numLevelSC);
                    tempIndex = findMortonIndex(tempMorton3, 0);
                }

                if (variables->tempParticlesPerDomain[tempIndex][1] >= variables->tempParticlesPerDomain[tempIndex][0]) {
                    variables->reallocateTempParticles(tempIndex, numParticles - i);
                }
                variables->tempParticles[tempIndex][variables->tempParticlesPerDomain[tempIndex][1]] = i + firstParticle;
                variables->tempParticlesPerDomain[tempIndex][1]++;
            }
        }
        else{
            // If not in bounds, place in dummy domain
            if (variables->tempParticlesPerDomain[numDomains][1] >= variables->tempParticlesPerDomain[numDomains][0]) {
                variables->reallocateTempParticles(numDomains, numParticles - i);
            }
            variables->tempParticles[numDomains][variables->tempParticlesPerDomain[numDomains][1]] = i + firstParticle;
            variables->tempParticlesPerDomain[numDomains][1]++;
        }
    }

    tempIndex = firstParticle;
    for (i = 0; i < numDomains; i++) {
        variables->domainLinkage[0][i].firstParticle = tempIndex;
        tempIndex += variables->tempParticlesPerDomain[i][1];
        variables->domainLinkage[0][i].lastParticle = tempIndex - 1;
    }

//    printf("xdims = %d, ydims = %d, zdims = %d\n", boxDimVals[0], boxDimVals[1], boxDimVals[2]);
//    printf("xdiff = %.9g, ydiff = %.9g, zdiff = %.9g\n", xdiff, ydiff, zdiff);
//    if (numLevelSC > 2) {
//        printf("test\n");
//    }
//    if (numDomains == 3) {
//        printf("test\n");
//    }
//        if (numDomains > 1){
//            printf("test\n");
//        }
    //    if (numDomains > 4){
    //        printf("test\n");
    //    }
}

void MultipoleMethodForceCalculation::setOptimumLevelTranslations(){
    //    float level_switch[7]={100,2e4,1.7e5,1.3e6,1e7,7e7,5e8}; // cpu-tree
    //    float level_switch[8]={10,100,1.3e4,1e5,7e5,5e6,3e7,1.5e8}; // cpu-fmm
    //    float level_switch[7]={minAvgPPB,1.3e4,1e5,7e5,5e6,3e7,1.5e8}; // cpu-fmm
    //    float level_switch[7]={1000,1e5,5e5,5e6,3e7,2e8,1.5e9}; // gpu-tree
    //    float level_switch[7] = {1000,1e5,7e5,7e6,5e7,3e8,2e9}; // gpu-fmm

    float level_switch[8] = {static_cast<float>(2 * FMMConstant::minAvgPPB), 1000, 1.3e4, 1e5, 7e5, 5e6, 3e7, 1.5e8};

    int domainMaxLevel, domainNumParticles;
    numBoxIndexFull = 0;
    maxLevel = 0;
    for (int i = 0; i < numDomains; i++) {
        domainMaxLevel = 0;
        domainNumParticles = variables->domainLinkage[0][i].lastParticle - variables->domainLinkage[0][i].firstParticle + 1;
        if (domainNumParticles < level_switch[0]) {
            domainMaxLevel = 0;
        } else if (domainNumParticles < level_switch[1]) {
            domainMaxLevel += 1;
        } else if (domainNumParticles < level_switch[2]) {
            domainMaxLevel += 2;
        } else if (domainNumParticles < level_switch[3]) {
            domainMaxLevel += 3;
        } else if (domainNumParticles < level_switch[4]) {
            domainMaxLevel += 4;
        } else if (domainNumParticles < level_switch[5]) {
            domainMaxLevel += 5;
        } else if (domainNumParticles < level_switch[6]) {
            domainMaxLevel += 6;
        } else if (domainNumParticles < level_switch[7]) {
            domainMaxLevel += 7;
        } else {
            domainMaxLevel += 8;
        }
        maxLevel = std::max(maxLevel, domainMaxLevel);
    }
    
    
    if (forcedMaxLevel != -2){
        maxLevel = forcedMaxLevel;
    }
    

    domainNumBoxIndexFull = pow(3, (3 * maxLevel));
    numBoxIndexFull = numDomains * domainNumBoxIndexFull;
//    printf("level   : %d\n", maxLevel);
}

void MultipoleMethodForceCalculation::presortMortonTranslations(){
    int i, ii, jj, nx, ny, nz, boxIndex;
    float boxSize;
    vec3<float> domainBoxMin;
    int domainNumParticles;
    vec3<int> tempBoxIndex3D;

//    boxSize = rootBoxSize / (1 << maxLevel);
    
    float numTotalBoxes = pow(3.0, maxLevel);
    
    boxSize = rootBoxSize / numTotalBoxes;

    for (i = 0; i < numDomains; i++) {
        MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][i].mortonIndex, tempBoxIndex3D);
        tempBoxIndex3D.x -= variables->offsetSC[0].x;
        tempBoxIndex3D.y -= variables->offsetSC[0].y;
        tempBoxIndex3D.z -= variables->offsetSC[0].z;

        domainNumParticles = variables->domainLinkage[0][i].lastParticle - variables->domainLinkage[0][i].firstParticle + 1;

        domainBoxMin.x = boxMin.x + tempBoxIndex3D.x * rootBoxSize;
        domainBoxMin.y = boxMin.y + tempBoxIndex3D.y * rootBoxSize;
        domainBoxMin.z = boxMin.z + tempBoxIndex3D.z * rootBoxSize;

        for (ii = 0; ii < domainNumParticles; ii++) {
            nx = int((particleSystem->chargedParticles.position[variables->tempParticles[i][ii]].x - domainBoxMin.x) / boxSize);
            ny = int((particleSystem->chargedParticles.position[variables->tempParticles[i][ii]].y - domainBoxMin.y) / boxSize);
            nz = int((particleSystem->chargedParticles.position[variables->tempParticles[i][ii]].z - domainBoxMin.z) / boxSize);

            //            if (nx >= (1 << maxLevel) || ny >= (1 << maxLevel) || nz >= (1 << maxLevel)){
            //                printf("test\n");
            //            }

            if (nx >= numTotalBoxes) nx--;
            if (ny >= numTotalBoxes) ny--;
            if (nz >= numTotalBoxes) nz--;

            boxIndex = 0;
            for (jj = 0; jj < maxLevel; jj++) {
                boxIndex += (nx % 3) * pow(3, (3 * jj + 1));
                nx /= 3;

                boxIndex += (ny % 3) * pow(3, (3 * jj));
                ny /= 3;

                boxIndex += (nz % 3) * pow(3, (3 * jj + 2));
                nz /= 3;
            }
            variables->domainMortonIndex[variables->domainLinkage[0][i].firstParticle + ii - firstParticle] = boxIndex;
        }
    }

    int outlierDomainNum = 0;
    int outlierOffset = lastParticle - numOutliers;
    for (i = 0; i < numOutliers; i++) {
        if (i == 0 || variables->outlierList[i].signedMorton != variables->outlierList[i - 1].signedMorton) {
            MultipoleMethodFunction::signedUnmorton3(variables->outlierList[i].signedMorton, tempBoxIndex3D);
//            if (variables->outlierList[i].signedMorton >= 27){
//                printf("test\n");
//            }
            tempBoxIndex3D.x -= variables->offsetSC[0].x;
            tempBoxIndex3D.y -= variables->offsetSC[0].y;
            tempBoxIndex3D.z -= variables->offsetSC[0].z;

            domainBoxMin.x = boxMin.x + tempBoxIndex3D.x * rootBoxSize;
            domainBoxMin.y = boxMin.y + tempBoxIndex3D.y * rootBoxSize;
            domainBoxMin.z = boxMin.z + tempBoxIndex3D.z * rootBoxSize;

            variables->outlierDomainLinkage[outlierDomainNum].firstParticle = outlierOffset + i;
            if (outlierDomainNum > 0) {
                variables->outlierDomainLinkage[outlierDomainNum - 1].lastParticle = outlierOffset + i - 1;
            }
            variables->outlierDomainLinkage[outlierDomainNum].mortonIndex = variables->outlierList[i].signedMorton;
            outlierDomainNum++;
        }

        nx = int((particleSystem->chargedParticles.position[variables->outlierList[i].index].x - domainBoxMin.x) / boxSize);
        ny = int((particleSystem->chargedParticles.position[variables->outlierList[i].index].y - domainBoxMin.y) / boxSize);
        nz = int((particleSystem->chargedParticles.position[variables->outlierList[i].index].z - domainBoxMin.z) / boxSize);

        boxIndex = 0;
        for (jj = 0; jj < maxLevel; jj++) {
            boxIndex += (nx % 3) * pow(3, (3 * jj + 1));
            nx /= 3;

            boxIndex += (ny % 3) * pow(3, (3 * jj));
            ny /= 3;

            boxIndex += (nz % 3) * pow(3, (3 * jj + 2));
            nz /= 3;
        }

        variables->domainMortonIndex[outlierOffset + i] = boxIndex;
    }
    if (outlierDomainNum > 0) {
        variables->outlierDomainLinkage[outlierDomainNum - 1].lastParticle = lastParticle - 1;
    }
}

void MultipoleMethodForceCalculation::mortonTranslations(){
    int i, ii, jj, nx, ny, nz, boxIndex, tempIndex;
    float boxSize;
    vec3<float> domainBoxMin;
    int domainNumParticles;
    vec3<int> tempBoxIndex3D;

//    boxSize = rootBoxSize / (1 << maxLevel);
    
    float numTotalBoxes = pow(3.0, maxLevel);
    
    boxSize = rootBoxSize / numTotalBoxes;

    for (i = 0; i < numDomains; i++) {
        MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][i].mortonIndex, tempBoxIndex3D);
        tempBoxIndex3D.x -= variables->offsetSC[0].x;
        tempBoxIndex3D.y -= variables->offsetSC[0].y;
        tempBoxIndex3D.z -= variables->offsetSC[0].z;

        domainNumParticles = variables->domainLinkage[0][i].lastParticle - variables->domainLinkage[0][i].firstParticle + 1;

        domainBoxMin.x = boxMin.x + tempBoxIndex3D.x * rootBoxSize;
        domainBoxMin.y = boxMin.y + tempBoxIndex3D.y * rootBoxSize;
        domainBoxMin.z = boxMin.z + tempBoxIndex3D.z * rootBoxSize;

        for (ii = variables->domainLinkage[0][i].firstParticle; ii <= variables->domainLinkage[0][i].lastParticle; ii++) {
            nx = int((particleSystem->chargedParticles.position[ii].x - domainBoxMin.x) / boxSize);
            ny = int((particleSystem->chargedParticles.position[ii].y - domainBoxMin.y) / boxSize);
            nz = int((particleSystem->chargedParticles.position[ii].z - domainBoxMin.z) / boxSize);

            if (nx >= numTotalBoxes) nx--;
            if (ny >= numTotalBoxes) ny--;
            if (nz >= numTotalBoxes) nz--;

            boxIndex = 0;
            for (jj = 0; jj < maxLevel; jj++) {
                boxIndex += (nx % 3) * pow(3, (3 * jj + 1));
                nx /= 3;

                boxIndex += (ny % 3) * pow(3, (3 * jj));
                ny /= 3;

                boxIndex += (nz % 3) * pow(3, (3 * jj + 2));
                nz /= 3;
            }
            variables->domainMortonIndex[ii - firstParticle] = boxIndex;
        }
    }

    int outlierOffset = lastParticle - numOutliers;
    for (i = 0; i < numOutlierDomains; i++) {
        for (ii = variables->outlierDomainLinkage[i].firstParticle; ii <= variables->outlierDomainLinkage[i].lastParticle; ii++) {
//            if (ii == 370){
//                printf("test\n");
//            }
//            
            if (ii == variables->outlierDomainLinkage[i].firstParticle) {
                MultipoleMethodFunction::signedUnmorton3(variables->outlierList[ii - outlierOffset].signedMorton, tempBoxIndex3D);
                tempBoxIndex3D.x -= variables->offsetSC[0].x;
                tempBoxIndex3D.y -= variables->offsetSC[0].y;
                tempBoxIndex3D.z -= variables->offsetSC[0].z;

                domainBoxMin.x = boxMin.x + tempBoxIndex3D.x * rootBoxSize;
                domainBoxMin.y = boxMin.y + tempBoxIndex3D.y * rootBoxSize;
                domainBoxMin.z = boxMin.z + tempBoxIndex3D.z * rootBoxSize;
            }
            nx = int((particleSystem->chargedParticles.position[ii].x - domainBoxMin.x) / boxSize);
            ny = int((particleSystem->chargedParticles.position[ii].y - domainBoxMin.y) / boxSize);
            nz = int((particleSystem->chargedParticles.position[ii].z - domainBoxMin.z) / boxSize);

            boxIndex = 0;
            for (jj = 0; jj < maxLevel; jj++) {
                boxIndex += (nx % 3) * pow(3, (3 * jj + 1));
                nx /= 3;

                boxIndex += (ny % 3) * pow(3, (3 * jj));
                ny /= 3;

                boxIndex += (nz % 3) * pow(3, (3 * jj + 2));
                nz /= 3;
            }

            variables->domainMortonIndex[ii - firstParticle] = boxIndex;
        }
    }
}

void MultipoleMethodForceCalculation::getBoxDataTranslations(int& numBoxIndex, 
        int toggle){
    int i, ii, currentIndex;
    int numParticles = lastParticle - firstParticle;

    if (toggle == 0) {
        mortonTranslations();
        numBoxIndex = 0;
        int domainBias;

        for (i = 0; i < numBoxIndexFull; i++) {
            variables->domainBoxIndexMask[i] = -1;
        }
        for (i = 0; i < numDomains; i++) {
            domainBias = i * domainNumBoxIndexFull;
            //        domainNumBoxIndex[i] = 0;
            currentIndex = -1;
            variables->domainOffset[0][i][maxLevel] = numBoxIndex;
            for (ii = variables->domainLinkage[0][i].firstParticle; ii <= variables->domainLinkage[0][i].lastParticle; ii++) {
                //            if (ii == 376){
                //                printf("test\n");
                //            }
//                if (ii == 303) {
//                    printf("test\n");
//                }
                if (variables->domainMortonIndex[ii - firstParticle] != currentIndex) {
                    variables->domainBoxIndexMask[domainBias + variables->domainMortonIndex[ii - firstParticle]] = numBoxIndex;
                    variables->domainBoxIndexFull[numBoxIndex] = variables->domainMortonIndex[ii - firstParticle];
                    variables->particleOffset[0][numBoxIndex] = ii;
                    if (numBoxIndex > 0) variables->particleOffset[1][numBoxIndex - 1] = ii - 1;
                    currentIndex = variables->domainMortonIndex[ii - firstParticle];
                    //                domainNumBoxIndex[i]++;
                    numBoxIndex++;
                }
            }
            variables->domainOffset[1][i][maxLevel] = numBoxIndex - 1;
        }
        variables->particleOffset[1][numBoxIndex - 1] = variables->domainLinkage[0][numDomains - 1].lastParticle;

        
        outlierNumBoxIndex = 0;
        for (i = 0; i < numOutlierDomains; i++) {
            currentIndex = -1;
            variables->outlierDomainOffset[0][i][maxLevel] = outlierNumBoxIndex;
            for (ii = variables->outlierDomainLinkage[i].firstParticle; ii <= variables->outlierDomainLinkage[i].lastParticle; ii++) {
                if (variables->domainMortonIndex[ii - firstParticle] != currentIndex) {
                    variables->outlierBoxIndexFull[outlierNumBoxIndex] = variables->domainMortonIndex[ii - firstParticle];
                    variables->outlierParticleOffset[0][outlierNumBoxIndex] = ii;
                    if (outlierNumBoxIndex > 0) variables->outlierParticleOffset[1][outlierNumBoxIndex - 1] = ii - 1;
                    currentIndex = variables->domainMortonIndex[ii - firstParticle];
                    outlierNumBoxIndex++;
                }
            }
            variables->outlierDomainOffset[1][i][maxLevel] = outlierNumBoxIndex - 1;
        }
        if (outlierNumBoxIndex > 0){
            variables->outlierParticleOffset[1][outlierNumBoxIndex - 1] = variables->outlierDomainLinkage[numOutlierDomains - 1].lastParticle;
        }
    } else if (toggle == 1) {
        numBoxIndex = 0;
        for (i = 0; i < numDomains; i++) {
            variables->domainBoxIndexMask[i] = -1;
            if (variables->domainLinkage[0][i].lastParticle >= variables->domainLinkage[0][i].firstParticle) {
                variables->domainBoxIndexMask[i] = numBoxIndex;
                variables->domainBoxIndexFull[numBoxIndex] = i;
                numBoxIndex++;
            }
        }

        outlierNumBoxIndex = numOutlierDomains;
    }
}

void MultipoleMethodForceCalculation::getBoxDataOfParentTranslations(int& numBoxIndex,
        int numLevel, int toggle){
    int i, ii, jj, numBoxIndexOld, currentIndex, boxIndex;

    int domainBias;

    if (toggle == 0) {

        if (numLevel == 0) variables->levelOffset[maxLevel] = variables->levelOffset[0] + numBoxIndex;
        else variables->levelOffset[numLevel - 1] = variables->levelOffset[numLevel] + numBoxIndex;

        for (i = 0; i < numBoxIndexFull; i++) {
            variables->domainBoxIndexMask[i] = -1;
        }
        
        float numLevelBoxes = pow(3, 3 * numLevel);

        numBoxIndexOld = numBoxIndex;
        numBoxIndex = 0;
        jj = 0;
        for (i = 0; i < numDomains; i++) {
            domainBias = i * numLevelBoxes;
            currentIndex = -1;
            variables->domainOffset[0][i][numLevel] = numBoxIndex;
            for (ii = variables->domainOffset[0][i][numLevel + 1]; ii <= variables->domainOffset[1][i][numLevel + 1]; ii++) {
                boxIndex = ii + variables->levelOffset[numLevel];
                if (currentIndex != variables->domainBoxIndexFull[boxIndex] / 27) {
                    currentIndex = variables->domainBoxIndexFull[boxIndex] / 27;
                    variables->domainBoxIndexMask[domainBias + currentIndex] = numBoxIndex;

                    if (numLevel == 0) {
                        variables->domainBoxIndexFull[numBoxIndex + variables->levelOffset[maxLevel]] = currentIndex;
                    } else {
                        variables->domainBoxIndexFull[numBoxIndex + variables->levelOffset[numLevel - 1]] = currentIndex;
                    }

                    if (treeOrFMM == 0) {
                        variables->particleOffset[0][numBoxIndex] = variables->particleOffset[0][jj];
                        if (numBoxIndex > 0) variables->particleOffset[1][numBoxIndex - 1] = variables->particleOffset[0][jj] - 1;
                    }
                    numBoxIndex++;
                }
                jj++;
            }
            variables->domainOffset[1][i][numLevel] = numBoxIndex - 1;
        }
        if (treeOrFMM == 0) variables->particleOffset[1][numBoxIndex - 1] = variables->particleOffset[1][numBoxIndexOld - 1];


        numBoxIndexOld = outlierNumBoxIndex;
        outlierNumBoxIndex = 0;
        jj = 0;
        for (i = 0; i < numOutlierDomains; i++) {
            currentIndex = -1;
            variables->outlierDomainOffset[0][i][numLevel] = outlierNumBoxIndex;
            for (ii = variables->outlierDomainOffset[0][i][numLevel + 1]; ii <= variables->outlierDomainOffset[1][i][numLevel + 1]; ii++) {
                if (currentIndex != variables->outlierBoxIndexFull[ii] / 27) {
                    currentIndex = variables->outlierBoxIndexFull[ii] / 27;
                    variables->outlierBoxIndexFull[outlierNumBoxIndex] = currentIndex;
                    variables->outlierParticleOffset[0][outlierNumBoxIndex] = variables->outlierParticleOffset[0][jj];
                    if (outlierNumBoxIndex > 0) variables->outlierParticleOffset[1][outlierNumBoxIndex - 1] = variables->outlierParticleOffset[0][jj] - 1;
                    outlierNumBoxIndex++;
                }
                jj++;
            }
            variables->outlierDomainOffset[1][i][numLevel] = outlierNumBoxIndex - 1;
        }
        variables->outlierParticleOffset[1][outlierNumBoxIndex - 1] = variables->outlierParticleOffset[1][numBoxIndexOld - 1];
    } else if (toggle == 1) {

        if (numLevel == 0) variables->levelOffset[maxLevel] = variables->levelOffset[0] + numBoxIndex;
        else variables->levelOffset[maxLevel + numLevel] = variables->levelOffset[numLevel + maxLevel - 1] + numBoxIndex;

        numBoxIndexOld = numBoxIndex;
        numBoxIndex = 0;
        for (i = 0; i < numBoxIndexFull; i++) {
            variables->domainBoxIndexMask[i] = -1;
        }

        for (i = 0; i < variables->numDomainsLevel[numLevel]; i++) {
            if (variables->domainLinkage[numLevel][i].lastParticle >= variables->domainLinkage[numLevel][i].firstParticle) {
                variables->domainBoxIndexMask[i] = numBoxIndex;
                variables->domainBoxIndexFull[numBoxIndex + variables->levelOffset[numLevel + maxLevel]] = i;
                numBoxIndex++;
            }
        }
        
        outlierNumBoxIndex = numOutlierDomains;
    }
}

void MultipoleMethodForceCalculation::getBoxIndexMaskTranslations(int numBoxIndex,
        int numLevel, int toggle){
    int i, ii, boxIndex;

    int domainBias;
    if (toggle == 0) {
        float numLevelBoxes = pow(3, (3 * numLevel));
        for (i = 0; i < numBoxIndexFull; i++) variables->domainBoxIndexMask[i] = -1;
        for (i = 0; i < numDomains; i++) {
            domainBias = i * numLevelBoxes;
            for (ii = variables->domainOffset[0][i][numLevel]; ii <= variables->domainOffset[1][i][numLevel]; ii++) {
                boxIndex = ii + variables->levelOffset[numLevel - 1];
                variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[boxIndex]] = ii;
            }
        }
    } else if (toggle == 1) {
        boxIndex = 0;
        for (i = 0; i < numBoxIndexFull; i++) variables->domainBoxIndexMask[i] = -1;
        for (i = 0; i < variables->numDomainsLevel[numLevel]; i++) {
            if (variables->domainLinkage[numLevel][i].lastParticle >= variables->domainLinkage[numLevel][i].firstParticle) {
                variables->domainBoxIndexMask[i] = boxIndex;
                boxIndex++;
            }
        }
    }
}

void MultipoleMethodForceCalculation::sortParticlesTranslations(){
    int i, ii, buff, domainNumParticles;
    int numParticles = lastParticle - firstParticle;

    presortMortonTranslations();

    // Sort particleSystem->chargedParticles within each domain using local morton index
    buff = 0;
    for (i = 0; i < numDomains; i++) {
        domainNumParticles = variables->domainLinkage[0][i].lastParticle - variables->domainLinkage[0][i].firstParticle + 1;

        for (ii = 0; ii < domainNumParticles; ii++) {
            //variables->sortValue[ii] = variables->domainMortonIndex[ii + buff];
            variables->sortValue[ii] = variables->domainMortonIndex[variables->domainLinkage[0][i].firstParticle + ii - firstParticle];
            variables->sortIndex[ii] = variables->tempParticles[i][ii];
        }
        sort(i, 0);
        for (ii = 0; ii < domainNumParticles; ii++) {
            variables->permutation[ii + buff] = variables->sortIndex[ii] - firstParticle;
        }
        buff += domainNumParticles;
    }

    int offset;
    for (i = 0; i < numOutlierDomains; i++) {
        offset = variables->outlierDomainLinkage[i].firstParticle - lastParticle + numOutliers;
        domainNumParticles = variables->outlierDomainLinkage[i].lastParticle - variables->outlierDomainLinkage[i].firstParticle + 1;
        for (ii = 0; ii < domainNumParticles; ii++) {
            variables->sortValue[ii] = variables->domainMortonIndex[variables->outlierDomainLinkage[i].firstParticle + ii - firstParticle];
            variables->sortIndex[ii] = variables->outlierList[ii + offset].index;
        }
        sortOutliers(i, 0);
        for (ii = 0; ii < domainNumParticles; ii++) {
            variables->permutation[ii + buff] = variables->sortIndex[ii] - firstParticle;
        }
        buff += domainNumParticles;
    }
    
    for (ii = 0; ii < variables->tempParticlesPerDomain[numDomains][1]; ii++){
        variables->permutation[ii + buff] = variables->tempParticles[numDomains][ii] - firstParticle;
    }

    vec4<float>* sortBuffer;
    sortBuffer = new vec4<float>[numParticles];
    for (i = 0; i < numParticles; i++) {
        sortBuffer[i] = particleSystem->chargedParticles.position[variables->permutation[i] + firstParticle];
    }
    for (i = 0; i < numParticles; i++) {
        particleSystem->chargedParticles.position[i + firstParticle] = sortBuffer[i];
    }
    for (i = 0; i < numParticles; i++) {
        sortBuffer[i] = particleSystem->chargedParticles.velocity[variables->permutation[i] + firstParticle];
    }
    for (i = 0; i < numParticles; i++) {
        particleSystem->chargedParticles.velocity[i + firstParticle] = sortBuffer[i];
    }
    for (i = 0; i < numParticles; i++) {
        sortBuffer[i] = particleSystem->chargedParticles.acceleration[variables->permutation[i] + firstParticle];
    }
    for (i = 0; i < numParticles; i++) {
        particleSystem->chargedParticles.acceleration[i + firstParticle] = sortBuffer[i];
    }
    
    delete[] sortBuffer;
    
    float* sortBuffer2 = new float[numParticles];
    for (i = 0; i < numParticles; i++) {
        sortBuffer2[i] = particleSystem->chargedParticles.temperature[variables->permutation[i] + firstParticle];
    }
    for (i = 0; i < numParticles; i++) {
        particleSystem->chargedParticles.temperature[i + firstParticle] = sortBuffer2[i];
    }
    for (i = 0; i < numParticles; i++){
        sortBuffer2[i] = particleSystem->chargedParticles.timeFlag[variables->permutation[i] + firstParticle];
    }
    for (i = 0; i < numParticles; i++){
        particleSystem->chargedParticles.timeFlag[i + firstParticle] = sortBuffer2[i];
    }
    
    delete[] sortBuffer2;
    
    int* sortBuffer3 = new int[numParticles];
    for (i = 0; i < numParticles; i++) {
        sortBuffer3[i] = particleSystem->chargedParticles.timesFragmented[variables->permutation[i] + firstParticle];
    }
    for (i = 0; i < numParticles; i++) {
        particleSystem->chargedParticles.timesFragmented[i + firstParticle] = sortBuffer3[i];
    }
    for (i = 0; i < numParticles; i++){
        sortBuffer3[i] = particleSystem->chargedParticles.crash[variables->permutation[i] + firstParticle];
    }
    for (i = 0; i < numParticles; i++){
        particleSystem->chargedParticles.crash[i + firstParticle] = sortBuffer3[i];
    }
    
    delete[] sortBuffer3;
}

void MultipoleMethodForceCalculation::countNonEmptyBoxesTranslations(){
    int i, ii, currentIndex, numLevel, domainNumParticles;
    int numParticles = lastParticle - firstParticle;

    mortonTranslations();

    // Sort particles
    int* tempSortValue = new int[numParticles];
    for (i = 0; i < numDomains; i++) {
        domainNumParticles = variables->domainLinkage[0][i].lastParticle - variables->domainLinkage[0][i].firstParticle + 1;
        for (ii = 0; ii < domainNumParticles; ii++) {
            variables->sortValue[ii] = variables->domainMortonIndex[variables->domainLinkage[0][i].firstParticle + ii - firstParticle];
            //            if (variables->domainLinkage[0][i].firstParticle + ii == 376 && numDomains == 5){
            //                printf("%d\n",variables->sortValue[ii]);
            //            }
            //            if (variables->domainLinkage[0][i].firstParticle + ii == 1049){
            //                printf("%d\n",variables->sortValue[ii]);
            //            }
            variables->sortIndex[ii] = variables->domainLinkage[0][i].firstParticle + ii;
        }
        sort(i, 1);
        for (ii = 0; ii < domainNumParticles; ii++) {
            tempSortValue[variables->domainLinkage[0][i].firstParticle + ii - firstParticle] = variables->sortValue[ii];
        }
    }

    // At each level, find how many boxes are filled using particles within each box
    numBoxIndexLeaf = 0;
    numBoxIndexTotal = 0;
    int checkIndex;
    for (numLevel = maxLevel; numLevel >= 1; numLevel--) {
        for (i = 0; i < numDomains; i++) {
            currentIndex = -1;
            domainNumParticles = variables->domainLinkage[0][i].lastParticle - variables->domainLinkage[0][i].firstParticle + 1;
            for (ii = 0; ii < domainNumParticles; ii++) {
                checkIndex = tempSortValue[variables->domainLinkage[0][i].firstParticle + ii - firstParticle] / (pow(3, 3 * (maxLevel - numLevel)));
                if (checkIndex != currentIndex) {
                    if (numLevel == maxLevel) numBoxIndexLeaf++;
                    numBoxIndexTotal++;
                    currentIndex = checkIndex;
                }
            }
        }
    }

    // Check for filled 0-th level supercells
    for (i = 0; i < numDomains; i++) {
        if (variables->domainLinkage[0][i].lastParticle >= variables->domainLinkage[0][i].firstParticle) {
            if (maxLevel == 0) numBoxIndexLeaf++;
            numBoxIndexTotal++;
        }
    }

    // Check for filled higher level supercells
    int tempBoxIndex, guessIndex, j, jj;
    int prevSCIndex, tempNumBoxIndex;
    vec3<int> tempBoxIndex3D;
    for (numLevel = 1; numLevel < numLevelSC; numLevel++) {
        currentIndex = -1;
        prevSCIndex = -1;
        tempNumBoxIndex = -1;
        for (i = 0; i < variables->numDomainsLevel[numLevel - 1]; i++) {
            if (currentIndex != variables->domainLinkage[numLevel - 1][i].mortonIndex / 27) {
                tempNumBoxIndex++;
                currentIndex = variables->domainLinkage[numLevel - 1][i].mortonIndex / 27;
                MultipoleMethodFunction::unmorton3(currentIndex, tempBoxIndex3D);
                tempBoxIndex3D.x += variables->offsetSC[numLevel].x;
                tempBoxIndex3D.y += variables->offsetSC[numLevel].y;
                tempBoxIndex3D.z += variables->offsetSC[numLevel].z;
                MultipoleMethodFunction::morton3(tempBoxIndex3D, tempBoxIndex, numLevelSC);

                // Find index of parent cell
                guessIndex = findMortonIndex(tempBoxIndex, numLevel, prevSCIndex);

                // Set first particle of domain equal to first particle in first child cell
                variables->domainLinkage[numLevel][guessIndex].firstParticle = variables->domainLinkage[numLevel - 1][i].firstParticle;

                // For every empty supercell, set first and last particles
                for (jj = prevSCIndex + 1; jj < guessIndex; jj++) {
                    variables->domainLinkage[numLevel][jj].firstParticle = variables->domainLinkage[numLevel - 1][i - 1].lastParticle;
                    variables->domainLinkage[numLevel][jj].lastParticle = variables->domainLinkage[numLevel - 1][i - 1].lastParticle - 1;
                }

                // Set last particle of previous domain
                if (i > 0) variables->domainLinkage[numLevel][prevSCIndex].lastParticle = variables->domainLinkage[numLevel - 1][i - 1].lastParticle;

                prevSCIndex = guessIndex;
                numBoxIndexTotal++;
            }
            // Assign parent cell for each supercell
            variables->domainLinkage[numLevel - 1][i].parent = tempNumBoxIndex;
        }
        variables->domainLinkage[numLevel][variables->numDomainsLevel[numLevel] - 1].lastParticle = variables->domainLinkage[numLevel - 1][variables->numDomainsLevel[numLevel - 1] - 1].lastParticle;
    }

    for (i = 0; i < numOutlierDomains; i++) {
        domainNumParticles = variables->outlierDomainLinkage[i].lastParticle - variables->outlierDomainLinkage[i].firstParticle + 1;
        for (ii = 0; ii < domainNumParticles; ii++) {
            variables->sortValue[ii] = variables->domainMortonIndex[variables->outlierDomainLinkage[i].firstParticle + ii - firstParticle];
            //            if (variables->domainLinkage[0][i].firstParticle + ii == 376 && numDomains == 5){
            //                printf("%d\n",variables->sortValue[ii]);
            //            }
            //            if (variables->domainLinkage[0][i].firstParticle + ii == 1049){
            //                printf("%d\n",variables->sortValue[ii]);
            //            }
            variables->sortIndex[ii] = variables->outlierDomainLinkage[i].firstParticle + ii;
        }
        sortOutliers(i, 1);
        for (ii = 0; ii < domainNumParticles; ii++) {
            tempSortValue[variables->outlierDomainLinkage[i].firstParticle + ii - firstParticle] = variables->sortValue[ii];
        }
    }

    int outlierNumBoxIndexLeaf = 0;
    for (i = 0; i < numOutlierDomains; i++) {
        currentIndex = -1;
        domainNumParticles = variables->outlierDomainLinkage[i].lastParticle - variables->outlierDomainLinkage[i].firstParticle + 1;
        for (ii = 0; ii < domainNumParticles; ii++) {
            checkIndex = tempSortValue[variables->outlierDomainLinkage[i].firstParticle + ii - firstParticle];
            if (checkIndex != currentIndex) {
                outlierNumBoxIndexLeaf++;
                currentIndex = checkIndex;
            }
        }
    }

    variables->outlierBoxIndexFull = new int[outlierNumBoxIndexLeaf];

    // Needed to ensure that interaction list can handle both regular cells and imaginary outlier cells
    numBoxIndexLeaf = std::max(numBoxIndexLeaf, outlierNumBoxIndexLeaf);

    delete[] tempSortValue;
}

void MultipoleMethodForceCalculation::getInteractionListTranslations(int numBoxIndex,
        int numLevel, FmmInteractionType interactionType){
    int jxmin, jxmax, jymin, jymax, jzmin, jzmax, ii, ib, jj, jb, ix, iy, iz, jx, jy, jz, boxIndex;
    int ixp, iyp, izp, jxp, jyp, jzp;
    vec3<int> boxIndex3D, domainIndex3D;
    int i, j, jdiff, halfDiff, jtemp;
    
    int xSym, ySym, zSym;
    
    int tempDomainIndex, tempQuad;
    vec3<int> symmetryCheck, symmetryHalfWidth;
    
    if (interactionType == FmmInteractionType::M2L_Subcell || 
            interactionType == FmmInteractionType::M2L_SubcellOutliers ||
            interactionType == FmmInteractionType::M2L_SubcellTopLevel ||
            interactionType == FmmInteractionType::M2L_SubcellTopLevelOutliers ||
            interactionType == FmmInteractionType::P2P ||
            interactionType == FmmInteractionType::P2P_Outliers){
        jdiff = pow(3, numLevel);
        symmetryHalfWidth.x = particleSystem->params.symmetryAxes.x * pow(3, particleSystem->params.numSymmetryLevels) / 2;
        symmetryHalfWidth.y = particleSystem->params.symmetryAxes.y * pow(3, particleSystem->params.numSymmetryLevels) / 2;
        symmetryHalfWidth.z = particleSystem->params.symmetryAxes.z * pow(3, particleSystem->params.numSymmetryLevels) / 2;
    }
    else if (interactionType == FmmInteractionType::M2L_Supercell ||
            interactionType == FmmInteractionType::M2L_SupercellOutliers ||
            interactionType == FmmInteractionType::M2L_SupercellTopLevel ||
            interactionType == FmmInteractionType::M2L_SupercellTopLevelOutliers ||
            interactionType == FmmInteractionType::P2P_Supercell ||
            interactionType == FmmInteractionType::P2P_SupercellOutliers){
        jdiff = pow(3, numLevel);
        symmetryHalfWidth.x = particleSystem->params.symmetryAxes.x * pow(3, particleSystem->params.numSymmetryLevels - numLevel) / 2;
        symmetryHalfWidth.y = particleSystem->params.symmetryAxes.y * pow(3, particleSystem->params.numSymmetryLevels - numLevel) / 2;
        symmetryHalfWidth.z = particleSystem->params.symmetryAxes.z * pow(3, particleSystem->params.numSymmetryLevels - numLevel) / 2;
    }
    halfDiff = (jdiff >> 1);

    int tempDomainNumBoxIndexFull = pow(3, (3 * numLevel));

    int domainBias, domainBoxIndex, guessIndex;
    vec3<int> tempDomainIndex3D, domainIndex3DHolder;

    // P2P
    if (interactionType == FmmInteractionType::P2P) {
        // Loop through every domain
        for (i = 0; i < numDomains; i++) {
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->domainOffset[0][i][numLevel]; ii <= variables->domainOffset[1][i][numLevel]; ii++) {
                ib = ii + variables->levelOffset[numLevel - 1];
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton3(variables->domainBoxIndexFull[ib], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;
                
                for (jx = ix - 1; jx <= ix + 1; jx++){
                    for (jy = iy - 1; jy <= iy + 1; jy++){
                        for (jz = iz - 1; jz <= iz + 1; jz++){
                            boxIndex3D.x = jx;
                            boxIndex3D.y = jy;
                            boxIndex3D.z = jz;
                            tempDomainIndex3D = domainIndex3D;
                            
                            symmetryCheck = 0;
                            
                            if (jx < 0) {
                                do {
                                    boxIndex3D.x += jdiff;
                                    tempDomainIndex3D.x--;
                                } while (boxIndex3D.x < 0);
                            } else if (jx >= jdiff) {
                                do {
                                    boxIndex3D.x -= jdiff;
                                    tempDomainIndex3D.x++;
                                } while (boxIndex3D.x >= jdiff);
                            }
                            
                            if (jy < 0) {
                                do {
                                    boxIndex3D.y += jdiff;
                                    tempDomainIndex3D.y--;
                                } while (boxIndex3D.y < 0);
                            } else if (jy >= jdiff) {
                                do {
                                    boxIndex3D.y -= jdiff;
                                    tempDomainIndex3D.y++;
                                } while (boxIndex3D.y >= jdiff);
                            }
                            
                            if (jz < 0) {
                                do {
                                    boxIndex3D.z += jdiff;
                                    tempDomainIndex3D.z--;
                                } while (boxIndex3D.z < 0);
                            } else if (jz >= jdiff) {
                                do {
                                    boxIndex3D.z -= jdiff;
                                    tempDomainIndex3D.z++;
                                } while (boxIndex3D.z >= jdiff);
                            }
                            
                            while (tempDomainIndex3D.x >= boxDimVals[0] + variables->offsetSC[0].x){
                                tempDomainIndex3D.x -= boxDimVals[0];
                                symmetryCheck.x++;
                            }
                            while (tempDomainIndex3D.x < variables->offsetSC[0].x){
                                tempDomainIndex3D.x += boxDimVals[0];
                                symmetryCheck.x--;
                            }
                            
                            while (tempDomainIndex3D.y >= boxDimVals[1] + variables->offsetSC[0].y){
                                tempDomainIndex3D.y -= boxDimVals[1];
                                symmetryCheck.y++;
                            }
                            while (tempDomainIndex3D.y < variables->offsetSC[0].y){
                                tempDomainIndex3D.y += boxDimVals[1];
                                symmetryCheck.y--;
                            }
                            
                            while (tempDomainIndex3D.z >= boxDimVals[2] + variables->offsetSC[0].z){
                                tempDomainIndex3D.z -= boxDimVals[2];
                                symmetryCheck.z++;
                            }
                            while (tempDomainIndex3D.z < variables->offsetSC[0].z){
                                tempDomainIndex3D.z += boxDimVals[2];
                                symmetryCheck.z--;
                            }
                            
                            // Check to see if domain is valid and cell contains any particles
                            if (abs(symmetryCheck.x) <= symmetryHalfWidth.x && abs(symmetryCheck.y) <= symmetryHalfWidth.y &&
                                    abs(symmetryCheck.z) <= symmetryHalfWidth.z) {
                                MultipoleMethodFunction::morton3(boxIndex3D, boxIndex, numLevel);
                                if (tempDomainIndex3D.x == domainIndex3D.x && tempDomainIndex3D.y == domainIndex3D.y && tempDomainIndex3D.z == domainIndex3D.z) {
                                    domainBias = i * tempDomainNumBoxIndexFull;
                                } else {
                                    MultipoleMethodFunction::morton3(tempDomainIndex3D, domainBoxIndex, numLevelSC);
                                    guessIndex = findMortonIndex(domainBoxIndex, 0);
                                    domainBias = guessIndex * tempDomainNumBoxIndexFull;
                                }
                                jj = variables->domainBoxIndexMask[boxIndex + domainBias];
                                if (jj != -1) {
                                    MultipoleMethodFunction::balancedTernary(symmetryCheck, tempDomainIndex);
                                    
                                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
                                    variables->interactionList[ii][variables->numInteraction[ii]][1] = tempDomainIndex;
                                    variables->numInteraction[ii]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }// P2P at domain level
    else if (interactionType == FmmInteractionType::P2P_Supercell) {
        // Loop through every domain
        //for (ii = 0; ii < numDomains; ii++){
        for (ii = 0; ii < numBoxIndex; ii++) {
            ib = ii + variables->levelOffset[numLevel + maxLevel];
            variables->numInteraction[ii] = 0;
            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][ii].mortonIndex, boxIndex3D);
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][variables->domainBoxIndexFull[ib]].mortonIndex, boxIndex3D);
            ix = boxIndex3D.x;
            iy = boxIndex3D.y;
            iz = boxIndex3D.z;
            // Check cells within a 3x3x3 box centered around box being considered
            for (jx = ix - 1; jx <= ix + 1; jx++) {
                for (jy = iy - 1; jy <= iy + 1; jy++) {
                    for (jz = iz - 1; jz <= iz + 1; jz++) {
                        symmetryCheck = 0;
                        
                        boxIndex3D.x = jx;
                        boxIndex3D.y = jy;
                        boxIndex3D.z = jz;
                        
                        while (boxIndex3D.x >= variables->dimValsSC[numLevel].x + variables->offsetSC[numLevel].x){
                            boxIndex3D.x -= variables->dimValsSC[numLevel].x;
                            symmetryCheck.x++;
                        }
                        while (boxIndex3D.x < variables->offsetSC[numLevel].x){
                            boxIndex3D.x += variables->dimValsSC[numLevel].x;
                            symmetryCheck.x--;
                        }

                        while (boxIndex3D.y >= variables->dimValsSC[numLevel].y + variables->offsetSC[numLevel].y){
                            boxIndex3D.y -= variables->dimValsSC[numLevel].y;
                            symmetryCheck.y++;
                        }
                        while (boxIndex3D.y < variables->offsetSC[numLevel].y){
                            boxIndex3D.y += variables->dimValsSC[numLevel].y;
                            symmetryCheck.y--;
                        }

                        while (boxIndex3D.z >= variables->dimValsSC[numLevel].z + variables->offsetSC[numLevel].z){
                            boxIndex3D.z -= variables->dimValsSC[numLevel].z;
                            symmetryCheck.z++;
                        }
                        while (boxIndex3D.z < variables->offsetSC[numLevel].z){
                            boxIndex3D.z += variables->dimValsSC[numLevel].z;
                            symmetryCheck.z--;
                        }
                        
                        if (abs(symmetryCheck.x) <= symmetryHalfWidth.x && abs(symmetryCheck.y) <= symmetryHalfWidth.y &&
                                abs(symmetryCheck.z) <= symmetryHalfWidth.z) {
                            MultipoleMethodFunction::morton3(boxIndex3D, domainBoxIndex, numLevelSC - numLevel);
                            guessIndex = findMortonIndex(domainBoxIndex, numLevel);
                            guessIndex = variables->domainBoxIndexMask[guessIndex];
                            if (guessIndex != -1) {
                                MultipoleMethodFunction::balancedTernary(symmetryCheck, tempDomainIndex);

                                variables->interactionList[ii][variables->numInteraction[ii]][0] = guessIndex;
                                variables->interactionList[ii][variables->numInteraction[ii]][1] = tempDomainIndex;
                                variables->numInteraction[ii]++;
                            }
                        }
                    }
                }
            }
        }
    }// P2P for outliers
    else if (interactionType == FmmInteractionType::P2P_Outliers) {
        // Loop through every domain
        for (i = 0; i < numOutlierDomains; i++) {
            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->outlierDomainOffset[0][i][numLevel]; ii <= variables->outlierDomainOffset[1][i][numLevel]; ii++) {
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton3(variables->outlierBoxIndexFull[ii], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;

                // Check cells within a 3x3x3 box centered around box being considered
                for (jx = ix - 1; jx <= ix + 1; jx++) {
                    for (jy = iy - 1; jy <= iy + 1; jy++) {
                        for (jz = iz - 1; jz <= iz + 1; jz++) {
                            boxIndex3D.x = jx;
                            boxIndex3D.y = jy;
                            boxIndex3D.z = jz;
                            tempDomainIndex3D = domainIndex3D;
                            
                            // Shift cell to lie within some domain
                            symmetryCheck = 0;
                            
                            if (jx < 0) {
                                do {
                                    boxIndex3D.x += jdiff;
                                    tempDomainIndex3D.x--;
                                } while (boxIndex3D.x < 0);
                            } else if (jx >= jdiff) {
                                do {
                                    boxIndex3D.x -= jdiff;
                                    tempDomainIndex3D.x++;
                                } while (boxIndex3D.x >= jdiff);
                            }
                            
                            if (jy < 0) {
                                do {
                                    boxIndex3D.y += jdiff;
                                    tempDomainIndex3D.y--;
                                } while (boxIndex3D.y < 0);
                            } else if (jy >= jdiff) {
                                do {
                                    boxIndex3D.y -= jdiff;
                                    tempDomainIndex3D.y++;
                                } while (boxIndex3D.y >= jdiff);
                            }
                            
                            if (jz < 0) {
                                do {
                                    boxIndex3D.z += jdiff;
                                    tempDomainIndex3D.z--;
                                } while (boxIndex3D.z < 0);
                            } else if (jz >= jdiff) {
                                do {
                                    boxIndex3D.z -= jdiff;
                                    tempDomainIndex3D.z++;
                                } while (boxIndex3D.z >= jdiff);
                            }
                            
                            while (tempDomainIndex3D.x >= boxDimVals[0] + variables->offsetSC[0].x){
                                tempDomainIndex3D.x -= boxDimVals[0];
                                symmetryCheck.x++;
                            }
                            while (tempDomainIndex3D.x < variables->offsetSC[0].x){
                                tempDomainIndex3D.x += boxDimVals[0];
                                symmetryCheck.x--;
                            }
                            
                            while (tempDomainIndex3D.y >= boxDimVals[1] + variables->offsetSC[0].y){
                                tempDomainIndex3D.y -= boxDimVals[1];
                                symmetryCheck.y++;
                            }
                            while (tempDomainIndex3D.y < variables->offsetSC[0].y){
                                tempDomainIndex3D.y += boxDimVals[1];
                                symmetryCheck.y--;
                            }
                            
                            while (tempDomainIndex3D.z >= boxDimVals[2] + variables->offsetSC[0].z){
                                tempDomainIndex3D.z -= boxDimVals[2];
                                symmetryCheck.z++;
                            }
                            while (tempDomainIndex3D.z < variables->offsetSC[0].z){
                                tempDomainIndex3D.z += boxDimVals[2];
                                symmetryCheck.z--;
                            }
                            
                            // Check to see if domain is valid and cell contains any particles
                            if (abs(symmetryCheck.x) <= symmetryHalfWidth.x && abs(symmetryCheck.y) <= symmetryHalfWidth.y &&
                                    abs(symmetryCheck.z) <= symmetryHalfWidth.z) {
                                MultipoleMethodFunction::morton3(boxIndex3D, boxIndex, numLevel);
                                if (tempDomainIndex3D.x == domainIndex3D.x && tempDomainIndex3D.y == domainIndex3D.y && tempDomainIndex3D.z == domainIndex3D.z) {
                                    domainBias = i * tempDomainNumBoxIndexFull;
                                } else {
                                    MultipoleMethodFunction::morton3(tempDomainIndex3D, domainBoxIndex, numLevelSC);
                                    guessIndex = findMortonIndex(domainBoxIndex, 0);
                                    domainBias = guessIndex * tempDomainNumBoxIndexFull;
                                }
                                jj = variables->domainBoxIndexMask[boxIndex + domainBias];
                                if (jj != -1) {
                                    MultipoleMethodFunction::balancedTernary(symmetryCheck, tempDomainIndex);
                                    
                                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
                                    variables->interactionList[ii][variables->numInteraction[ii]][1] = tempDomainIndex;
                                    variables->numInteraction[ii]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }// P2P at domain level for outliers
    else if (interactionType == FmmInteractionType::P2P_SupercellOutliers) {
        // Loop through every domain
        //for (ii = 0; ii < numDomains; ii++){
        for (ii = 0; ii < numBoxIndex; ii++) {
            variables->numInteraction[ii] = 0;
            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][ii].mortonIndex, boxIndex3D);
            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[ii].mortonIndex, boxIndex3D);
            ix = boxIndex3D.x;
            iy = boxIndex3D.y;
            iz = boxIndex3D.z;
            // Check cells within a 3x3x3 box centered around box being considered
            for (jx = ix - 1; jx <= ix + 1; jx++) {
                for (jy = iy - 1; jy <= iy + 1; jy++) {
                    for (jz = iz - 1; jz <= iz + 1; jz++) {
                        symmetryCheck = 0;
                        
                        boxIndex3D.x = jx;
                        boxIndex3D.y = jy;
                        boxIndex3D.z = jz;
                        
                       while (boxIndex3D.x >= variables->dimValsSC[numLevel].x + variables->offsetSC[numLevel].x){
                            boxIndex3D.x -= variables->dimValsSC[numLevel].x;
                            symmetryCheck.x++;
                        }
                        while (boxIndex3D.x < variables->offsetSC[numLevel].x){
                            boxIndex3D.x += variables->dimValsSC[numLevel].x;
                            symmetryCheck.x--;
                        }

                        while (boxIndex3D.y >= variables->dimValsSC[numLevel].y + variables->offsetSC[numLevel].y){
                            boxIndex3D.y -= variables->dimValsSC[numLevel].y;
                            symmetryCheck.y++;
                        }
                        while (boxIndex3D.y < variables->offsetSC[numLevel].y){
                            boxIndex3D.y += variables->dimValsSC[numLevel].y;
                            symmetryCheck.y--;
                        }

                        while (boxIndex3D.z >= variables->dimValsSC[numLevel].z + variables->offsetSC[numLevel].z){
                            boxIndex3D.z -= variables->dimValsSC[numLevel].z;
                            symmetryCheck.z++;
                        }
                        while (boxIndex3D.z < variables->offsetSC[numLevel].z){
                            boxIndex3D.z += variables->dimValsSC[numLevel].z;
                            symmetryCheck.z--;
                        }
                        
                        if (abs(symmetryCheck.x) <= symmetryHalfWidth.x && abs(symmetryCheck.y) <= symmetryHalfWidth.y &&
                                abs(symmetryCheck.z) <= symmetryHalfWidth.z) {
                            MultipoleMethodFunction::morton3(boxIndex3D, domainBoxIndex, numLevelSC - numLevel);
                            guessIndex = findMortonIndex(domainBoxIndex, numLevel);
                            guessIndex = variables->domainBoxIndexMask[guessIndex];
                            if (guessIndex != -1) {
                                MultipoleMethodFunction::balancedTernary(symmetryCheck, tempDomainIndex);

                                variables->interactionList[ii][variables->numInteraction[ii]][0] = guessIndex;
                                variables->interactionList[ii][variables->numInteraction[ii]][1] = tempDomainIndex;
                                variables->numInteraction[ii]++;
                            }
                        }
                    }
                }
            }
        }
    }// M2L at level 1
    else if (interactionType == FmmInteractionType::M2L_SubcellTopLevel) {
        // Loop through every domain
        for (i = 0; i < numDomains; i++) {
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->domainOffset[0][i][numLevel]; ii <= variables->domainOffset[1][i][numLevel]; ii++) {
                ib = ii + variables->levelOffset[numLevel - 1];
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton3(variables->domainBoxIndexFull[ib], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;

                // Loop through every domain
                for (j = 0; j < numDomains; j++) {
                    MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][j].mortonIndex, tempDomainIndex3D);
                    domainBias = j * pow(3, (3 * numLevel));

                    // Loop through all cells within domain j
                    for (jj = variables->domainOffset[0][j][numLevel]; jj <= variables->domainOffset[1][j][numLevel]; jj++) {
                        jb = jj + variables->levelOffset[numLevel - 1];
                        MultipoleMethodFunction::unmorton3(variables->domainBoxIndexFull[jb], boxIndex3D);
//                        jx = boxIndex3D.x;
//                        jy = boxIndex3D.y;
//                        jz = boxIndex3D.z;
//                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);  // Should add zero
//                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);  // Should add zero
//                        jz += jdiff * (tempDomainIndex3D.z - domainIndex3D.z);

                        for (xSym = -(symmetryHalfWidth.x != 0); xSym <= (symmetryHalfWidth.x != 0); xSym++){
                            for (ySym = -(symmetryHalfWidth.y != 0); ySym <= (symmetryHalfWidth.y != 0); ySym++){
                                for (zSym = -(symmetryHalfWidth.z != 0); zSym <= (symmetryHalfWidth.z != 0); zSym++){
                                    domainIndex3DHolder = tempDomainIndex3D;
                                    domainIndex3DHolder.x += xSym;
                                    domainIndex3DHolder.y += ySym;
                                    domainIndex3DHolder.z += zSym;
                                    jx = boxIndex3D.x;
                                    jy = boxIndex3D.y;
                                    jz = boxIndex3D.z;
                                    jx += jdiff * (domainIndex3DHolder.x - domainIndex3D.x);
                                    jy += jdiff * (domainIndex3DHolder.y - domainIndex3D.y);
                                    jz += jdiff * (domainIndex3DHolder.z - domainIndex3D.z);
                                    
                                    symmetryCheck.x = xSym;
                                    symmetryCheck.y = ySym;
                                    symmetryCheck.z = zSym;
                                    
                                    if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                                        MultipoleMethodFunction::balancedTernary(symmetryCheck, tempDomainIndex);
                                        variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
                                        variables->interactionList[ii][variables->numInteraction[ii]][1] = tempDomainIndex;
                                        variables->numInteraction[ii]++;
                                    }
                                }
                            }
                        }
//                        // Check if jj-th cell lies within 3x3x3 box of ii-th cell
//                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
//                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
//                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 0;
//                            variables->numInteraction[ii]++;
//                        }
                        
                    }
                }
            }
        }
    }// M2L at level 1 for outliers
    else if (interactionType == FmmInteractionType::M2L_SubcellTopLevelOutliers) {
//        // Loop through every domain
//        for (i = 0; i < numOutlierDomains; i++) {
//            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[i].mortonIndex, domainIndex3D);
//            // Loop through all cells within domain i
//            for (ii = variables->outlierDomainOffset[0][i][numLevel]; ii <= variables->outlierDomainOffset[1][i][numLevel]; ii++) {
//                variables->numInteraction[ii] = 0;
//
//                // Find local box index
//                MultipoleMethodFunction::unmorton(variables->outlierBoxIndexFull[ii], boxIndex3D);
//                ix = boxIndex3D.x;
//                iy = boxIndex3D.y;
//                iz = boxIndex3D.z;
//
//                // Loop through every domain
//                for (j = 0; j < numDomains; j++) {
//                    MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][j].mortonIndex, tempDomainIndex3D);
//                    domainBias = j * (1 << (3 * numLevel));
//
//                    // Loop through all cells within domain j
//                    for (jj = variables->domainOffset[0][j][numLevel]; jj <= variables->domainOffset[1][j][numLevel]; jj++) {
//                        jb = jj + variables->levelOffset[numLevel - 1];
//                        MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[jb], boxIndex3D);
//                        jx = boxIndex3D.x;
//                        jy = boxIndex3D.y;
//                        jz = boxIndex3D.z;
//                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
//                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
//                        jz += jdiff * (tempDomainIndex3D.z - domainIndex3D.z);
//
//                        // Check if jj-th cell lies within 3x3x3 box of ii-th cell
//                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
//                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
//                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 0;
//                            variables->numInteraction[ii]++;
//                        }
//                        
//                        // Second quadrant
//                        jy = boxIndex3D.x;
//                        jx = jdiff - boxIndex3D.y - 1;
//                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
//                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
//                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
//                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
//                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 1;
//                            variables->numInteraction[ii]++;
//                        }
//                        
//                        // Third quadrant
//                        jx = jdiff - boxIndex3D.x - 1;
//                        jy = jdiff - boxIndex3D.y - 1;
//                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
//                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
//                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
//                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
//                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 2;
//                            variables->numInteraction[ii]++;
//                        }
//                        
//                        // Fourth quadrant
//                        jx = boxIndex3D.y;
//                        jy = jdiff - boxIndex3D.x - 1;
//                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
//                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
//                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
//                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
//                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 3;
//                            variables->numInteraction[ii]++;
//                        }
//                    }
//                }
//            }
//        }
    }// M2L at subcell level
    else if (interactionType == FmmInteractionType::M2L_Subcell) {
        // Loop through every domain
        for (i = 0; i < numDomains; i++) {
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->domainOffset[0][i][numLevel]; ii <= variables->domainOffset[1][i][numLevel]; ii++) {
                ib = ii + variables->levelOffset[numLevel - 1];
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton3(variables->domainBoxIndexFull[ib], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;
                ixp = (ix + 3) / 3;
                iyp = (iy + 3) / 3;
                izp = (iz + 3) / 3;

                // Loop through all cells within a 3x3x3 box centered around the parent cell
                for (jxp = ixp - 1; jxp <= ixp + 1; jxp++) {
                    for (jyp = iyp - 1; jyp <= iyp + 1; jyp++) {
                        for (jzp = izp - 1; jzp <= izp + 1; jzp++) {
                            for (jx = 3 * jxp - 3; jx <= 3 * jxp - 1; jx++) {
                                for (jy = 3 * jyp - 3; jy <= 3 * jyp - 1; jy++) {
                                    for (jz = 3 * jzp - 3; jz <= 3 * jzp - 1; jz++) {
                                        // Cells must be outside of a 3x3x3 box centered around ii-th cell
                                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                                            boxIndex3D.x = jx;
                                            boxIndex3D.y = jy;
                                            boxIndex3D.z = jz;
                                            tempDomainIndex3D = domainIndex3D;
                                            
                                            // Shift cell to lie within some domain
                                            symmetryCheck = 0;

                                            if (jx < 0) {
                                                do {
                                                    boxIndex3D.x += jdiff;
                                                    tempDomainIndex3D.x--;
                                                } while (boxIndex3D.x < 0);
                                            } else if (jx >= jdiff) {
                                                do {
                                                    boxIndex3D.x -= jdiff;
                                                    tempDomainIndex3D.x++;
                                                } while (boxIndex3D.x >= jdiff);
                                            }

                                            if (jy < 0) {
                                                do {
                                                    boxIndex3D.y += jdiff;
                                                    tempDomainIndex3D.y--;
                                                } while (boxIndex3D.y < 0);
                                            } else if (jy >= jdiff) {
                                                do {
                                                    boxIndex3D.y -= jdiff;
                                                    tempDomainIndex3D.y++;
                                                } while (boxIndex3D.y >= jdiff);
                                            }

                                            if (jz < 0) {
                                                do {
                                                    boxIndex3D.z += jdiff;
                                                    tempDomainIndex3D.z--;
                                                } while (boxIndex3D.z < 0);
                                            } else if (jz >= jdiff) {
                                                do {
                                                    boxIndex3D.z -= jdiff;
                                                    tempDomainIndex3D.z++;
                                                } while (boxIndex3D.z >= jdiff);
                                            }

                                            while (tempDomainIndex3D.x >= boxDimVals[0] + variables->offsetSC[0].x){
                                                tempDomainIndex3D.x -= boxDimVals[0];
                                                symmetryCheck.x++;
                                            }
                                            while (tempDomainIndex3D.x < variables->offsetSC[0].x){
                                                tempDomainIndex3D.x += boxDimVals[0];
                                                symmetryCheck.x--;
                                            }

                                            while (tempDomainIndex3D.y >= boxDimVals[1] + variables->offsetSC[0].y){
                                                tempDomainIndex3D.y -= boxDimVals[1];
                                                symmetryCheck.y++;
                                            }
                                            while (tempDomainIndex3D.y < variables->offsetSC[0].y){
                                                tempDomainIndex3D.y += boxDimVals[1];
                                                symmetryCheck.y--;
                                            }

                                            while (tempDomainIndex3D.z >= boxDimVals[2] + variables->offsetSC[0].z){
                                                tempDomainIndex3D.z -= boxDimVals[2];
                                                symmetryCheck.z++;
                                            }
                                            while (tempDomainIndex3D.z < variables->offsetSC[0].z){
                                                tempDomainIndex3D.z += boxDimVals[2];
                                                symmetryCheck.z--;
                                            }

                                            if (abs(symmetryCheck.x) <= symmetryHalfWidth.x && abs(symmetryCheck.y) <= symmetryHalfWidth.y &&
                                                    abs(symmetryCheck.z) <= symmetryHalfWidth.z) {
                                                MultipoleMethodFunction::morton3(boxIndex3D, boxIndex, numLevel);
                                                if (tempDomainIndex3D.x == domainIndex3D.x && tempDomainIndex3D.y == domainIndex3D.y && tempDomainIndex3D.z == domainIndex3D.z) {
                                                    domainBias = i * tempDomainNumBoxIndexFull;
                                                } else {
                                                    MultipoleMethodFunction::morton3(tempDomainIndex3D, domainBoxIndex, numLevelSC);
                                                    guessIndex = findMortonIndex(domainBoxIndex, 0);
                                                    domainBias = guessIndex * tempDomainNumBoxIndexFull;
                                                }
                                                jj = variables->domainBoxIndexMask[boxIndex + domainBias];
                                                if (jj != -1) {
                                                    MultipoleMethodFunction::balancedTernary(symmetryCheck, tempDomainIndex);

                                                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
                                                    variables->interactionList[ii][variables->numInteraction[ii]][1] = tempDomainIndex;
                                                    variables->numInteraction[ii]++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }// M2L at subcell level for outliers
    else if (interactionType == FmmInteractionType::M2L_SubcellOutliers) {
        // Loop through every domain
        for (i = 0; i < numOutlierDomains; i++) {
            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->outlierDomainOffset[0][i][numLevel]; ii <= variables->outlierDomainOffset[1][i][numLevel]; ii++) {
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton3(variables->outlierBoxIndexFull[ii], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;
                ixp = (ix + 3) / 3;
                iyp = (iy + 3) / 3;
                izp = (iz + 3) / 3;

                // Loop through all cells within a 3x3x3 box centered around the parent cell
                for (jxp = ixp - 1; jxp <= ixp + 1; jxp++) {
                    for (jyp = iyp - 1; jyp <= iyp + 1; jyp++) {
                        for (jzp = izp - 1; jzp <= izp + 1; jzp++) {
                            for (jx = 3 * jxp - 3; jx <= 3 * jxp - 1; jx++) {
                                for (jy = 3 * jyp - 3; jy <= 3 * jyp - 1; jy++) {
                                    for (jz = 3 * jzp - 3; jz <= 3 * jzp - 1; jz++) {
                                        // Cells must be outside of a 3x3x3 box centered around ii-th cell
                                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                                            boxIndex3D.x = jx;
                                            boxIndex3D.y = jy;
                                            boxIndex3D.z = jz;
                                            tempDomainIndex3D = domainIndex3D;
                                            
                                            // Shift cell to lie within some domain
                                            symmetryCheck = 0;

                                            if (jx < 0) {
                                                do {
                                                    boxIndex3D.x += jdiff;
                                                    tempDomainIndex3D.x--;
                                                } while (boxIndex3D.x < 0);
                                            } else if (jx >= jdiff) {
                                                do {
                                                    boxIndex3D.x -= jdiff;
                                                    tempDomainIndex3D.x++;
                                                } while (boxIndex3D.x >= jdiff);
                                            }

                                            if (jy < 0) {
                                                do {
                                                    boxIndex3D.y += jdiff;
                                                    tempDomainIndex3D.y--;
                                                } while (boxIndex3D.y < 0);
                                            } else if (jy >= jdiff) {
                                                do {
                                                    boxIndex3D.y -= jdiff;
                                                    tempDomainIndex3D.y++;
                                                } while (boxIndex3D.y >= jdiff);
                                            }

                                            if (jz < 0) {
                                                do {
                                                    boxIndex3D.z += jdiff;
                                                    tempDomainIndex3D.z--;
                                                } while (boxIndex3D.z < 0);
                                            } else if (jz >= jdiff) {
                                                do {
                                                    boxIndex3D.z -= jdiff;
                                                    tempDomainIndex3D.z++;
                                                } while (boxIndex3D.z >= jdiff);
                                            }
                                            
                                            while (tempDomainIndex3D.x >= boxDimVals[0] + variables->offsetSC[0].x){
                                                tempDomainIndex3D.x -= boxDimVals[0];
                                                symmetryCheck.x++;
                                            }
                                            while (tempDomainIndex3D.x < variables->offsetSC[0].x){
                                                tempDomainIndex3D.x += boxDimVals[0];
                                                symmetryCheck.x--;
                                            }

                                            while (tempDomainIndex3D.y >= boxDimVals[1] + variables->offsetSC[0].y){
                                                tempDomainIndex3D.y -= boxDimVals[1];
                                                symmetryCheck.y++;
                                            }
                                            while (tempDomainIndex3D.y < variables->offsetSC[0].y){
                                                tempDomainIndex3D.y += boxDimVals[1];
                                                symmetryCheck.y--;
                                            }

                                            while (tempDomainIndex3D.z >= boxDimVals[2] + variables->offsetSC[0].z){
                                                tempDomainIndex3D.z -= boxDimVals[2];
                                                symmetryCheck.z++;
                                            }
                                            while (tempDomainIndex3D.z < variables->offsetSC[0].z){
                                                tempDomainIndex3D.z += boxDimVals[2];
                                                symmetryCheck.z--;
                                            }

                                            if (abs(symmetryCheck.x) <= symmetryHalfWidth.x && abs(symmetryCheck.y) <= symmetryHalfWidth.y &&
                                                    abs(symmetryCheck.z) <= symmetryHalfWidth.z) {
                                                MultipoleMethodFunction::morton3(boxIndex3D, boxIndex, numLevel);
                                                if (tempDomainIndex3D.x == domainIndex3D.x && tempDomainIndex3D.y == domainIndex3D.y && tempDomainIndex3D.z == domainIndex3D.z) {
                                                    domainBias = i * tempDomainNumBoxIndexFull;
                                                } else {
                                                    MultipoleMethodFunction::morton3(tempDomainIndex3D, domainBoxIndex, numLevelSC);
                                                    guessIndex = findMortonIndex(domainBoxIndex, 0);
                                                    domainBias = guessIndex * tempDomainNumBoxIndexFull;
                                                }
                                                jj = variables->domainBoxIndexMask[boxIndex + domainBias];
                                                if (jj != -1) {
                                                    MultipoleMethodFunction::balancedTernary(symmetryCheck, tempDomainIndex);

                                                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
                                                    variables->interactionList[ii][variables->numInteraction[ii]][1] = tempDomainIndex;
                                                    variables->numInteraction[ii]++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }// M2L for supercells
    else if (interactionType == FmmInteractionType::M2L_Supercell) {
        // Set bounds for supercell indices
        jxmin = variables->offsetSC[numLevel].x;
        jxmax = variables->offsetSC[numLevel].x + variables->dimValsSC[numLevel].x - 1;
        jymin = variables->offsetSC[numLevel].y;
        jymax = variables->offsetSC[numLevel].y + variables->dimValsSC[numLevel].y - 1;
        jzmin = variables->offsetSC[numLevel].z;
        jzmax = variables->offsetSC[numLevel].z + variables->dimValsSC[numLevel].z - 1;

        // Loop through every domain at level numLevel
        //for (ii = 0; ii < variables->numDomainsLevel[numLevel]; ii++){
        for (ii = 0; ii < numBoxIndex; ii++) {
            ib = ii + variables->levelOffset[maxLevel + numLevel];
            variables->numInteraction[ii] = 0;

            // Find box index of domain
            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][ii].mortonIndex, boxIndex3D);
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].mortonIndex, boxIndex3D);
            ix = boxIndex3D.x;
            iy = boxIndex3D.y;
            iz = boxIndex3D.z;
            ixp = (ix + 3) / 3;
            iyp = (iy + 3) / 3;
            izp = (iz + 3) / 3;

            // Loop through all cells within a 3x3x3 box centered around the parent cell
            for (jxp = ixp - 1; jxp <= ixp + 1; jxp++) {
                for (jyp = iyp - 1; jyp <= iyp + 1; jyp++) {
                    for (jzp = izp - 1; jzp <= izp + 1; jzp++) {
                        for (jx = 3 * jxp - 3; jx <= 3 * jxp - 1; jx++) {
                            for (jy = 3 * jyp - 3; jy <= 3 * jyp - 1; jy++) {
                                for (jz = 3 * jzp - 3; jz <= 3 * jzp - 1; jz++) {
                                    // Cells must be outside of a 3x3x3 box centered around ii-th cell
                                    if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                                        boxIndex3D.x = jx;
                                        boxIndex3D.y = jy;
                                        boxIndex3D.z = jz;
                                        
                                        symmetryCheck = 0;
                                        
                                        while (boxIndex3D.x >= variables->dimValsSC[numLevel].x + variables->offsetSC[numLevel].x){
                                            boxIndex3D.x -= variables->dimValsSC[numLevel].x;
                                            symmetryCheck.x++;
                                        }
                                        while (boxIndex3D.x < variables->offsetSC[numLevel].x){
                                            boxIndex3D.x += variables->dimValsSC[numLevel].x;
                                            symmetryCheck.x--;
                                        }

                                        while (boxIndex3D.y >= variables->dimValsSC[numLevel].y + variables->offsetSC[numLevel].y){
                                            boxIndex3D.y -= variables->dimValsSC[numLevel].y;
                                            symmetryCheck.y++;
                                        }
                                        while (boxIndex3D.y < variables->offsetSC[numLevel].y){
                                            boxIndex3D.y += variables->dimValsSC[numLevel].y;
                                            symmetryCheck.y--;
                                        }

                                        while (boxIndex3D.z >= variables->dimValsSC[numLevel].z + variables->offsetSC[numLevel].z){
                                            boxIndex3D.z -= variables->dimValsSC[numLevel].z;
                                            symmetryCheck.z++;
                                        }
                                        while (boxIndex3D.z < variables->offsetSC[numLevel].z){
                                            boxIndex3D.z += variables->dimValsSC[numLevel].z;
                                            symmetryCheck.z--;
                                        }

                                        if (abs(symmetryCheck.x) <= symmetryHalfWidth.x && abs(symmetryCheck.y) <= symmetryHalfWidth.y &&
                                                abs(symmetryCheck.z) <= symmetryHalfWidth.z) {
                                            MultipoleMethodFunction::morton3(boxIndex3D, domainBoxIndex, numLevelSC - numLevel);
                                            guessIndex = findMortonIndex(domainBoxIndex, numLevel);
                                            guessIndex = variables->domainBoxIndexMask[guessIndex];
                                            if (guessIndex != -1) {
                                                MultipoleMethodFunction::balancedTernary(symmetryCheck, tempDomainIndex);

                                                variables->interactionList[ii][variables->numInteraction[ii]][0] = guessIndex;
                                                variables->interactionList[ii][variables->numInteraction[ii]][1] = tempDomainIndex;
                                                variables->numInteraction[ii]++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }// M2L at supercell level for outliers
    else if (interactionType == FmmInteractionType::M2L_SupercellOutliers) {
        // Set bounds for supercell indices
        jxmin = variables->offsetSC[numLevel].x;
        jxmax = variables->offsetSC[numLevel].x + variables->dimValsSC[numLevel].x - 1;
        jymin = variables->offsetSC[numLevel].y;
        jymax = variables->offsetSC[numLevel].y + variables->dimValsSC[numLevel].y - 1;
        jzmin = variables->offsetSC[numLevel].z;
        jzmax = variables->offsetSC[numLevel].z + variables->dimValsSC[numLevel].z - 1;
        
        if (numOutliers == 3){
            printf("test\n");
        }

        // Loop through every domain at level numLevel
        //for (ii = 0; ii < variables->numDomainsLevel[numLevel]; ii++){
        for (ii = 0; ii < numBoxIndex; ii++) {
            variables->numInteraction[ii] = 0;

            // Find box index of domain
            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][ii].mortonIndex, boxIndex3D);
            //            MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].mortonIndex, boxIndex3D);
            //            ix = boxIndex3D.x;
            //            iy = boxIndex3D.y;
            //            iz = boxIndex3D.z;
            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[ii].mortonIndex, boxIndex3D);
            ix = boxIndex3D.x;
            iy = boxIndex3D.y;
            iz = boxIndex3D.z;

            if (ix < variables->boxMinIndexSC[numLevel].x && ((ix - variables->boxMinIndexSC[numLevel].x) % jdiff != 0)) {
                ix = (ix - variables->boxMinIndexSC[numLevel].x) / jdiff - 1;
            } else {
                ix = (ix - variables->boxMinIndexSC[numLevel].x) / jdiff;
            }

            if (iy < variables->boxMinIndexSC[numLevel].y && ((iy - variables->boxMinIndexSC[numLevel].y) % jdiff != 0)) {
                iy = (iy - variables->boxMinIndexSC[numLevel].y) / jdiff - 1;
            } else {
                iy = (iy - variables->boxMinIndexSC[numLevel].y) / jdiff;
            }

            if (iz < variables->boxMinIndexSC[numLevel].z && ((iz - variables->boxMinIndexSC[numLevel].z) % jdiff != 0)) {
                iz = (iz - variables->boxMinIndexSC[numLevel].z) / jdiff - 1;
            } else {
                iz = (iz - variables->boxMinIndexSC[numLevel].z) / jdiff;
            }

            ixp = (ix + 3) / 3;
            iyp = (iy + 3) / 3;
            izp = (iz + 3) / 3;

            // Loop through all cells within a 3x3x3 box centered around the parent cell
            for (jxp = ixp - 1; jxp <= ixp + 1; jxp++) {
                for (jyp = iyp - 1; jyp <= iyp + 1; jyp++) {
                    for (jzp = izp - 1; jzp <= izp + 1; jzp++) {
                        for (jx = 3 * jxp - 3; jx <= 3 * jxp - 1; jx++) {
                            for (jy = 3 * jyp - 3; jy <= 3 * jyp - 1; jy++) {
                                for (jz = 3 * jzp - 3; jz <= 3 * jzp - 1; jz++) {
                                    // Cells must be outside of a 3x3x3 box centered around ii-th cell
                                    if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                                        boxIndex3D.x = jx;
                                        boxIndex3D.y = jy;
                                        boxIndex3D.z = jz;
                                        
                                        symmetryCheck = 0;
                                        
                                        while (boxIndex3D.x >= variables->dimValsSC[numLevel].x + variables->offsetSC[numLevel].x){
                                            boxIndex3D.x -= variables->dimValsSC[numLevel].x;
                                            symmetryCheck.x++;
                                        }
                                        while (boxIndex3D.x < variables->offsetSC[numLevel].x){
                                            boxIndex3D.x += variables->dimValsSC[numLevel].x;
                                            symmetryCheck.x--;
                                        }

                                        while (boxIndex3D.y >= variables->dimValsSC[numLevel].y + variables->offsetSC[numLevel].y){
                                            boxIndex3D.y -= variables->dimValsSC[numLevel].y;
                                            symmetryCheck.y++;
                                        }
                                        while (boxIndex3D.y < variables->offsetSC[numLevel].y){
                                            boxIndex3D.y += variables->dimValsSC[numLevel].y;
                                            symmetryCheck.y--;
                                        }

                                        while (boxIndex3D.z >= variables->dimValsSC[numLevel].z + variables->offsetSC[numLevel].z){
                                            boxIndex3D.z -= variables->dimValsSC[numLevel].z;
                                            symmetryCheck.z++;
                                        }
                                        while (boxIndex3D.z < variables->offsetSC[numLevel].z){
                                            boxIndex3D.z += variables->dimValsSC[numLevel].z;
                                            symmetryCheck.z--;
                                        }

                                        if (abs(symmetryCheck.x) <= symmetryHalfWidth.x && abs(symmetryCheck.y) <= symmetryHalfWidth.y &&
                                                abs(symmetryCheck.z) <= symmetryHalfWidth.z) {
                                            MultipoleMethodFunction::morton3(boxIndex3D, domainBoxIndex, numLevelSC - numLevel);
                                            guessIndex = findMortonIndex(domainBoxIndex, numLevel);
                                            guessIndex = variables->domainBoxIndexMask[guessIndex];
                                            if (guessIndex != -1) {
                                                MultipoleMethodFunction::balancedTernary(symmetryCheck, tempDomainIndex);

                                                variables->interactionList[ii][variables->numInteraction[ii]][0] = guessIndex;
                                                variables->interactionList[ii][variables->numInteraction[ii]][1] = tempDomainIndex;
                                                variables->numInteraction[ii]++;
                                                if (numLevel == 1){
                                                    printf("test\n");
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }// M2L for supercells
    else if (interactionType == FmmInteractionType::M2L_SupercellTopLevel) {
        // Should be equivalent to regular M2L_Supercell. Main benefit in the
        // normal case is to only have to check domains that you know exist, but
        // with potential translational symmetry conditions, you'd have to check
        // all domains and all possible translational symmetry options, which
        // seems potentially worse than the regular supercell case.
        
        
//        // Loop through every domain at level numLevel
//        //for (ii = 0; ii < variables->numDomainsLevel[numLevel]; ii++){
//        for (ii = 0; ii < numBoxIndex; ii++) {
//            ib = ii + variables->levelOffset[maxLevel + numLevel];
//            variables->numInteraction[ii] = 0;
//
//            // Find box index of domain
//            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][ii].mortonIndex, boxIndex3D);
//            MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].mortonIndex, boxIndex3D);
//            ix = boxIndex3D.x;
//            iy = boxIndex3D.y;
//            iz = boxIndex3D.z;
//
//            // Loop through every domain at level numLevel
//            //for (jj = 0; jj < variables->numDomainsLevel[numLevel]; jj++){
//            for (jj = 0; jj < numBoxIndex; jj++) {
//                jb = jj + variables->levelOffset[maxLevel + numLevel];
//                //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][jj].mortonIndex, boxIndex3D);
//                MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex, boxIndex3D);
//                jx = boxIndex3D.x;
//                jy = boxIndex3D.y;
//                jz = boxIndex3D.z;
//                
//                
//                
//                if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
//                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
//                    variables->numInteraction[ii]++;
//                }
//            }
//        }
    }// M2L for supercells for outliers
    else if (interactionType == FmmInteractionType::M2L_SupercellTopLevelOutliers) {
//        // Loop through every domain at level numLevel
//
//        //        vec3<int> netOffset = {0, 0, 0};
//        //        int powThree = pow(3, numLevel);
//        //        for (ii = 0; ii < numLevel; ii++){
//        //            netOffset.x += powThree * variables->offsetSC[ii].x;
//        //            netOffset.y += powThree * variables->offsetSC[ii].y;
//        //            netOffset.z += powThree * variables->offsetSC[ii].z;
//        //            powThree *= 3;
//        //        }
//        //int jOffset = jdiff * 3 / 2;
//
//        //for (ii = 0; ii < variables->numDomainsLevel[numLevel]; ii++){
//        for (ii = 0; ii < outlierNumBoxIndex; ii++) {
//            variables->numInteraction[ii] = 0;
//
//            // Find box index of domain
//            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][ii].mortonIndex, boxIndex3D);
//            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[ii].mortonIndex, boxIndex3D);
//            ix = boxIndex3D.x;
//            iy = boxIndex3D.y;
//            iz = boxIndex3D.z;
//
//            if (ix < variables->boxMinIndexSC[numLevel].x) {
//                ix = (ix - variables->boxMinIndexSC[numLevel].x) / jdiff - 1;
//            } else {
//                ix = (ix - variables->boxMinIndexSC[numLevel].x) / jdiff;
//            }
//
//            if (iy < variables->boxMinIndexSC[numLevel].y) {
//                iy = (iy - variables->boxMinIndexSC[numLevel].y) / jdiff - 1;
//            } else {
//                iy = (iy - variables->boxMinIndexSC[numLevel].y) / jdiff;
//            }
//
//            if (iz < variables->boxMinIndexSC[numLevel].z) {
//                iz = (iz - variables->boxMinIndexSC[numLevel].z) / jdiff - 1;
//            } else {
//                iz = (iz - variables->boxMinIndexSC[numLevel].z) / jdiff;
//            }
//
//
//            // Loop through every domain at level numLevel
//            //for (jj = 0; jj < variables->numDomainsLevel[numLevel]; jj++){
//            for (jj = 0; jj < numBoxIndex; jj++) {
//                //                jb = jj + variables->levelOffset[maxLevel + numLevel];
//                //                //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][jj].mortonIndex, boxIndex3D);
//                //                MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex, boxIndex3D);
//                //                jx = variables->boxMinIndexSC[i].x + jdiff * boxIndex3D.x + jdiff / 2;
//                //                jy = variables->boxMinIndexSC[i].y + jdiff * boxIndex3D.y + jdiff / 2;
//                //                jz = variables->boxMinIndexSC[i].z + jdiff * boxIndex3D.z + jdiff / 2;
//                //                if (jx < ix - jOffset || ix + jOffset < jx || jy < iy - jOffset || iy + jOffset < jy || jz < iz - jOffset || iz + jOffset < jz){
//                //                    variables->interactionList[ii][variables->numInteraction[ii]] = jj;
//                //                    variables->numInteraction[ii]++;
//                //                }
//                jb = jj + variables->levelOffset[maxLevel + numLevel];
//                //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][jj].mortonIndex, boxIndex3D);
//                MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex, boxIndex3D);
//                jx = boxIndex3D.x;
//                jy = boxIndex3D.y;
//                jz = boxIndex3D.z;
//                if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
//                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
//                    variables->numInteraction[ii]++;
//                }
//            }
//            //            if (numDomains == 3){
//            //                printf("M2P variables->numInteraction[%d] = %d\n",ii,variables->numInteraction[ii]);
//            //            }
//        }
//    }
    }
}