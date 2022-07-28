/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "simulation_parameters.h"

//SimulationParameters::SimulationParameters(float timeStep, float simTime,
//        int numItersStored) {
//    timeStep = timeStep;
//    simTime = simTime;
//    numItersStored = numItersStored;
//    
//}

///**
// * Close save files and free info
// */
//SimulationParameters::~SimulationParameters(){
//    switch(saveDataToggle){
//        case (SaveDataOptions::original):
//            fclose(saveDataFiles[0]);
//            break;
//        case (SaveDataOptions::compressed):
//            fclose(saveDataFiles[0]);
//            break;
//        case (SaveDataOptions::compressedNew):
//            fclose(saveDataFiles[0]);
//            break;
//        case (SaveDataOptions::bothOld):
//            fclose(saveDataFiles[0]);
//            fclose(saveDataFiles[1]);
//            break;
//        case (SaveDataOptions::bothNew):
//            fclose(saveDataFiles[0]);
//            fclose(saveDataFiles[1]);
//            break;
//    }
//    delete[] saveDataFiles;
//}

///**
// * Initialize the electric field interpolation scheme
// * 
// * @param fieldsName File name for electric field data
// * @param trianglesName File name for triangulation points
// */
//void SimulationParameters::initializeInterpolationScheme(const char* fieldsName,
//        const char* trianglesName){
//    printf("Computing Delaunay triangulation (needed for interpolation)\n");
//    
//    interpScheme.initializeInterpolationScheme(fieldsName, trianglesName);
//}

/**
 * Deconstructor for SimulationParameters
 * Deallocate memory
 */
SimulationParameters::~SimulationParameters(){
    for (int i = 0; i < numSaveDataFiles; i++){
        fclose(saveDataFiles[i]);
    }
    delete[] saveDataFiles;
    
    if (eventsFile != nullptr) fclose(eventsFile);
    if (neutralsFile != nullptr) fclose(neutralsFile);
    
    freeInjectionTimes();
}

/**
 * Reads simulation conditions from external file
 * @param infoFileName Name of external file to read simulation info from
 */
void SimulationParameters::readInputInfoFile(const char* infoFileName){
    int i;
    
    char* infoName = concat(folder, infoFileName);
    
    float* info = new float[InfoConstant::INFO_SIZE];
    
    FILE* infop = fopen(infoName, "r");
    if (!infop){
        fileOpenError(infoName, __LINE__, __FILE__);
    }
    for (i = 0; i < InfoConstant::INFO_SIZE; i++){
        fscanf(infop, "%f", &info[i]);
    }
    fclose(infop);
    
    I0 = info[0];
    r0 = info[1];
    E0 = info[2];
    v0 = info[3];
    
    worldDimR = info[5];
    worldDimZ = info[6];
    
    a = info[7];
    nu0 = info[8];
    V = info[9];
    
    geom[0] = info[10]; // extractor distance
    geom[1] = info[11]; // extractor width
    geom[2] = info[12]; // extractor thickness
    geom[3] = info[13]; // outbound limit
    geom[4] = info[14]; // interpolation boundary
    
    geom[5] = info[16]; // offset    
    
    free(info);
    
    // Assume constantly injecting for entire simulation
    setDefaultInjectionTimes();
    
    setTimeStepInfo();
}

/**
 * Compute particles per timestep and total number of particles
 */
void SimulationParameters::setTimeStepInfo(){
    dN = timeStep / (PhysicalConstant::ELEM_CHARGE / I0);
    numIters = roundf(simTime / timeStep);
    interval = std::max(numIters / numItersStored, 1);
    
//    numParticles = ceil(dN * numIters);
    numParticles = ceil(dN * totalInjectingTime / timeStep);
}

/**
 * Create an output file of information to be used during post-processing
 */
void SimulationParameters::createOutputInfoFile(){
    if (I0 == 0.0){
        printf("Warning: Current not initialized!\n");
    }
    
    const float vz = InfoConstant::vz;
    float vAvg = vz + E0 * timeStep * PhysicalConstant::ELEM_CHARGE / 
            (PhysicalConstant::AMU * PhysicalConstant::EMI_MASS * 2);
    
    int InA = roundf(I0 * 1e9);
    
    FILE* fp;
    
    char infoNameTemp[512];
    snprintf(infoNameTemp, sizeof (char) * 512, "%dnA_%.0fps_information.txt", InA, 
            timeStep * 1e12); //creating filename depending on current

    char* infoNameStore = concat(folder, infoNameTemp);
    
    setTimeStepInfo();
    
    // Does nothing???
    int timeStepFactor = 10; // Factor between smallest and largest timeStep, specify here
    if ((timeStepFactor & 1) == 1) {
        printf("Error: Time step factor must be an even value!\n");
        exit(1);
    }
    
    fp = fopen(infoNameStore, "w+");
    fprintf(fp, "%i\n%i\n%i\n%i\n%i\n%f\n%i\n%i\n%f\n%f\n%f\n%f\n%f\n%f\n%i\n", InA, 
            numParticles, 3 * numParticles, numIters, numItersStored, timeStep * 1e9,
            InfoConstant::SAVE_SIZE_IONS, InfoConstant::SAVE_SIZE_NEUTRALS, r0, vAvg, worldDimZ, 
            geom[0], geom[1], geom[2], interval);
    // Current(nA), Number of particles, Size of particle array, Number of steps, Number of steps stored, Time interval (ns), Size of the body struct, Injection radius(nm), Injection velocity (m/s), Dimensional Voltage [V], extractor distance, extractor radius, extractor thickness, interval
    
    free(infoNameStore);
    fclose(fp);
}

/**
 * Set the injection start and stop times to continuously inject particles
 *      throughout the entire simulation
 */
void SimulationParameters::setDefaultInjectionTimes(){
    const int size = 1;
    freeInjectionTimes();
    startInjectionTimes = new float[size];
    stopInjectionTimes = new float[size];
    sizeInjectionTimes = size;
    
    startInjectionTimes[0] = 0.0;
    stopInjectionTimes[0] = simTime;
    
    totalInjectingTime = simTime;
}

/**
 * Set the injection start and stop times for the simulation. Can be used to 
 *      inject a group of particles, wait, inject more particles, etc.
 * @param startTimes - array of times where particles will start being injected
 * @param stopTimes - array of times where particles will stop being injected
 * @param size - size of the startTimes and stopTimes arrays
 */
void SimulationParameters::setInjectionTimes(float* startTimes, float* stopTimes, int size){
    int i;
    freeInjectionTimes();
    startInjectionTimes = new float[size];
    stopInjectionTimes = new float[size];
    sizeInjectionTimes = size;
    
    totalInjectingTime = 0;
    
    for (i = 0; i < size; i++){
        startInjectionTimes[i] = startTimes[i];
        stopInjectionTimes[i] = stopTimes[i];
        
        if (i > 0){
            if (startInjectionTimes[i] <= startInjectionTimes[i - 1]){
                printf("Error: Start injection times not strictly increasing.\n");
                exit(1);
            }
            if (stopInjectionTimes[i] <= stopInjectionTimes[i - 1]){
                printf("Error: Stop injection times not strictly increasing.\n");
                exit(1);
            }
        }
        if (stopInjectionTimes[i] <= startInjectionTimes[i]){
            printf("Error: Stop injection time must be greater than start injection time.\n");
            exit(1);
        }
        
        totalInjectingTime += (stopInjectionTimes[i] - startInjectionTimes[i]);
    }
}

/**
 * Initialize the state file based on the current, time step, and sim time
 * @param fileName - name of file to save state data to
 * @return - pointer to state data file
 */
FILE* SimulationParameters::initializeWritingStateFiles(const char* fileName){
    char bufferTemp[512];
    int InA = roundf(I0 / 1e-9);
    if (timeStep * 1e12 >= 0.999){
        snprintf(bufferTemp, sizeof(char) * 512, "state_%d_%.0fps_%.0fps", InA, timeStep*1e12, simTime * 1e12); // creating filename depending on current
    }
    else {
        snprintf(bufferTemp, sizeof(char) * 512, "state_%d_%.0ffs_%.0fps", InA, timeStep*1e15, simTime * 1e12); // creating filename depending on current
    }
    char* buffer2;
    char* buffer;
    buffer2 = concat(bufferTemp,fileName);
    buffer = concat(folder, buffer2);
    if (remove(buffer) == 0) {
            printf("\nOverwriting previous %s. Hope you saved the data!\n", fileName);
          }
          else {
            printf("\nCreating %s and storing your results.\n", fileName);
          }
    FILE* fileToWrite = fopen(buffer, "ab");
    free(buffer);
    free(buffer2);
    return fileToWrite;
}

/**
 * Initialize the event file based on the current, time step, and sim time
 * @param fileName - name of file to save event data to
 * @return - pointer to event data file
 */
FILE* SimulationParameters::initializeWritingEventsFiles(const char* fileName){
    char bufferTemp[512];
    int InA = roundf(I0 / 1e-9);
    if (timeStep * 1e12 >= 0.999){
        snprintf(bufferTemp, sizeof(char) * 512, "events_%d_%.0fps_%.0fps", InA, timeStep*1e12, simTime * 1e12); // creating filename depending on current
    }
    else {
        snprintf(bufferTemp, sizeof(char) * 512, "events_%d_%.0ffs_%.0fps", InA, timeStep*1e15, simTime * 1e12); // creating filename depending on current
    }
    char* buffer2;
    char* buffer;
    buffer2 = concat(bufferTemp,fileName);
    buffer = concat(folder, buffer2);
    if (remove(buffer) == 0) {
            printf("\nOverwriting previous %s. Hope you saved the data!\n", fileName);
          }
          else {
            printf("\nCreating %s and storing your results.\n", fileName);
          }
    FILE* fileToWrite = fopen(buffer, "ab");
    free(buffer);
    free(buffer2);
    return fileToWrite;
}

/**
 * Initialize the neutrals file based on the current, time step, and sim time
 * @param fileName - name of file to save neutrals data to
 * @return - pointer to neutrals data file
 */
FILE* SimulationParameters::initializeWritingNeutralsFiles(const char* fileName){
    char bufferTemp[512];
    int InA = roundf(I0 / 1e-9);
    if (timeStep * 1e12 >= 0.999){
        snprintf(bufferTemp, sizeof(char) * 512, "neutrals_%d_%.0fps_%.0fps", InA, timeStep*1e12, simTime * 1e12); // creating filename depending on current
    }
    else {
        snprintf(bufferTemp, sizeof(char) * 512, "neutrals_%d_%.0ffs_%.0fps", InA, timeStep*1e15, simTime * 1e12); // creating filename depending on current
    }
    char* buffer2;
    char* buffer;
    buffer2 = concat(bufferTemp,fileName);
    buffer = concat(folder, buffer2);
    if (remove(buffer) == 0) {
            printf("\nOverwriting previous %s. Hope you saved the data!\n", fileName);
          }
          else {
            printf("\nCreating %s and storing your results.\n", fileName);
          }
    FILE* fileToWrite = fopen(buffer, "ab");
    free(buffer);
    free(buffer2);
    return fileToWrite;
}

/**
 * Initialize the data storage files
 * @param additionalName - optional additional name to differentiate files from
 *      one another
 */
void SimulationParameters::initializeDataStorageFiles(const char* additionalName){
    if (saveDataFiles != nullptr){
        for (int i = 0; i < numSaveDataFiles; i++){
            fclose(saveDataFiles[i]);
        }
        delete[] saveDataFiles;
    }
    if (eventsFile != nullptr) fclose(eventsFile);
    if (neutralsFile != nullptr) fclose(neutralsFile);
    
    switch (saveDataToggle){
        case (SaveDataOptions::original):
            saveDataFiles = new FILE*[1];
            numSaveDataFiles = 1;
            if (additionalName[0] == '\0'){
                saveDataFiles[0] = initializeWritingStateFiles("_truncated.bin");
            }
            else {
                char* buffer = concat(additionalName, "_truncated.bin");
                saveDataFiles[0] = initializeWritingStateFiles(buffer);
            }
            break;
        case (SaveDataOptions::newOriginal):
            saveDataFiles = new FILE*[1];
            numSaveDataFiles = 1;
            if (additionalName[0] == '\0'){
                saveDataFiles[0] = initializeWritingStateFiles("_truncatedNew.bin");
            }
            else{
                char* buffer = concat(additionalName, "_truncatedNew.bin");
                saveDataFiles[0] = initializeWritingStateFiles(buffer);
            }
            break;
    }
    
    if (additionalName[0] == '\0'){
        eventsFile = initializeWritingEventsFiles(".bin");
        neutralsFile = initializeWritingNeutralsFiles(".bin");
    }
    else {
        char* buffer = concat(additionalName, ".bin");
        eventsFile = initializeWritingEventsFiles(buffer);
        neutralsFile = initializeWritingNeutralsFiles(buffer);
    }
}