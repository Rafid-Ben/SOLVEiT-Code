/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include "interpolation_scheme.h"

/**
 * Destructor for InterpolationScheme
 * Deallocate all memory in the interpolation
 */
InterpolationScheme::~InterpolationScheme(){
    freeInterpolationScheme();
};

/**
 * Deallocate all memory associated with interpolation scheme
 */
void InterpolationScheme::freeInterpolationScheme(){
    if (m_elementNeighbor != nullptr) delete[] m_elementNeighbor;
    if (m_triangle != nullptr) delete[] m_triangle;
    if (m_nodeXY != nullptr) delete[] m_nodeXY;
    if (m_zEr != nullptr) delete[] m_zEr;
    if (m_zEz != nullptr) delete[] m_zEz;
    
    m_elementNeighbor = nullptr;
    m_triangle = nullptr;
    m_nodeXY = nullptr;
    m_zEr = nullptr;
    m_zEz = nullptr;
}

/**
 * Set the reference position where the electric field is zero
 * 
 * Useful for energy calculation reference points
 */
void InterpolationScheme::setReferencePosition(){
    switch (bgField) {
        case (BackgroundEfield::none):
            referencePosition = 0.0f;
            break;
        case (BackgroundEfield::singleEmitter):
            referencePosition.z = 140 * r0;
            // Position where Efield reaches zero and becomes negative in interpolation
            break;
    }
}

/**
 * Set the background electric field type
 * @param bgField - background electric field type
 */
void InterpolationScheme::changeBackgroundField(BackgroundEfield bgField) {
    if (bgField != this->bgField){
        freeInterpolationScheme();
        this->bgField = bgField;
        setReferencePosition();
    }
}

/**
 * Fill the interpolation scheme to be used
 * @param eFieldName - name of file containing interpolation scheme data
 */
void InterpolationScheme::fillInterpolationScheme(char* eFieldName) {
    FILE* fp = fopen(eFieldName, "r");
    if (!fp) {
        fileOpenError(eFieldName, __LINE__, __FILE__);
    }

    char buf[1024];
    int fieldCount = 0;
    int i;

    for (i = 0; i < m_nodeNum; i++) {
        fgets(buf, 1024, fp);
        fieldCount = 0;
        char *field = strtok(buf, ",");
        while (field) {
            if (fieldCount == 0) {
                m_nodeXY[2 * i] = atof(field);
            }
            if (fieldCount == 1) {
                m_nodeXY[2 * i + 1] = atof(field);
            }
            if (fieldCount == 2) {
                m_zEr[i] = atof(field);
            }
            if (fieldCount == 3) {
                m_zEz[i] = atof(field);
            }

            field = strtok(nullptr, ",");
            fieldCount++;
        }
    }
    fclose(fp);
}

/**
 * Fill the triangulation scheme to be used
 * @param triangleName - name of file containing triangulation scheme data
 */
void InterpolationScheme::fillTriangleScheme(char* triangleName) {
    FILE* fp = fopen(triangleName, "r");
    if (!fp) {
        fileOpenError(triangleName, __LINE__, __FILE__);
    }

    char buf[1024];
    int fieldCount = 0;
    int i;

    for (i = 0; i < m_elementNum; i++) {
        fgets(buf, 1024, fp);
        fieldCount = 0;
        char *field = strtok(buf, ",");
        while (field) {
            if (fieldCount == 0) {
                m_triangle[3 * i] = atoi(field);
            }
            if (fieldCount == 1) {
                m_triangle[3 * i + 1] = atoi(field);
            }
            if (fieldCount == 2) {
                m_triangle[3 * i + 2] = atoi(field);
            }
            if (fieldCount == 3) {
                m_elementNeighbor[3 * i] = atoi(field);
            }
            if (fieldCount == 4) {
                m_elementNeighbor[3 * i + 1] = atoi(field);
            }
            if (fieldCount == 5) {
                m_elementNeighbor[3 * i + 2] = atoi(field);
            }

            field = strtok(nullptr, ",");
            fieldCount++;
        }
    }
    fclose(fp);
}

/**
 * Initialize the interpolation scheme for a single emitter
 * @param folderName - name of folder containing necessary interpolation data
 * @param fieldsName - name of file containing interpolation scheme data
 * @param trianglesName - name of file containing triangulation scheme data
 */
void InterpolationScheme::initializeSingleEmitterInterpolationScheme(const char* folderName,
        const char* fieldsName, const char* trianglesName) {
    
    printf("Computing Delaunay triangulation (needed for interpolation)\n");

    char* eFieldName = concat(folderName, fieldsName);
    char* triangleName = concat(folderName, trianglesName);

    m_elementNum = returnSize(triangleName);
    m_nodeNum = returnSize(eFieldName);

    freeInterpolationScheme();
    m_elementNeighbor = new int[3 * m_elementNum];
    m_triangle = new int[3 * m_elementNum];
    m_nodeXY = new double[2 * m_nodeNum];
    m_zEr = new double[m_nodeNum];
    m_zEz = new double[m_nodeNum];

    fillInterpolationScheme(eFieldName);
    fillTriangleScheme(triangleName);
}

/**
 * Initialize the desired interpolation scheme
 * @param folderName - name of folder containing necessary interpolation data
 * @param fieldsName - name of file containing interpolation scheme data
 * @param trianglesName - name of file containing triangulation scheme data
 */
void InterpolationScheme::initializeInterpolationScheme(const char* folderName, 
        const char* fieldsName, const char* trianglesName) {
    
    switch (this->bgField){
        case (BackgroundEfield::singleEmitter):
            initializeSingleEmitterInterpolationScheme(folderName, fieldsName, trianglesName);
            break;
        case (BackgroundEfield::none):
            printf("No background field selected. Continuing without interpolation.\n");
            break;
    }
}

/**
 * Compute the background electric field for the single emitter interpolation scheme
 * @param x - x position coordinate
 * @param y - y position coordinate
 * @param z - z position coordinate
 * @return - electric field vector at the given position
 */
vec3<float> InterpolationScheme::computeEfieldSingleEmitter(float x, float y, float z) {
    if (m_elementNeighbor == nullptr || m_triangle == nullptr || m_nodeXY == nullptr || 
            m_zEr == nullptr || m_zEz == nullptr){
        printf("\nError: Interpolation scheme for single emitter not initialized!\n");
        exit(1);
    }
    float phi;
    phi = std::atan(y / x);
    if (x < 0){
        phi += PhysicalConstant::PI;
    }
    if (x == 0 && y != 0){
        phi = PhysicalConstant::PI / 2;
    }
    else if (x == 0 && y == 0){
        phi = 0;
    }
    
    // Normalize position
    x /= r0;
    y /= r0;
    z /= r0;
    
    vec3<float> Efield;
    
//    double phi;
//    float phi;
    
    double *zi;
    double *ri;
    double coords[2];
    coords[0] = (double) sqrt((double)(x * x + y * y));
    coords[1] = (double) z;
    
//    coords[0] = 0.0032663240563124418;
    
    if (z < 0 || z > 140){
        // Out of interpolation bounds
        Efield.x = 0;
        Efield.y = 0;
        Efield.z = 0;
    }
    else {
        // Compute polar interpolation
        ri = pwl_interp_2d_scattered_value(m_nodeNum, m_nodeXY, m_zEr,
                m_elementNum, m_triangle, m_elementNeighbor, 1, coords);

        zi = pwl_interp_2d_scattered_value(m_nodeNum, m_nodeXY, m_zEz,
                m_elementNum, m_triangle, m_elementNeighbor, 1, coords);
        
//        phi = atan(y / x);
//        if (x < 0){
//            phi += PhysicalConstant::PI;
//        }
        
        Efield.x = (float) *ri;
        Efield.x *= std::cos((float)phi);
        
        Efield.y = (float) *ri;
        Efield.y *= std::sin((float)phi);
        
        Efield.z = (float) *zi;
        
        free(ri);
        free(zi);
    }
    
    // Scale Efield
    Efield *= E0;
    
    return Efield;
}

/**
 * Evaluate the electric field for the desired interpolation scheme
 * @param x - x position coordinate
 * @param y - y position coordinate
 * @param z - z position coordinate
 * @return - electric field vector at the given position
 */
vec3<float> InterpolationScheme::Efield(float x, float y, float z) {
    vec3<float> Efield = 0.0;
    switch (this->bgField){
        case (BackgroundEfield::singleEmitter):
            Efield = computeEfieldSingleEmitter(x, y, z);
            break;
        case (BackgroundEfield::none):
            Efield = 0.0;
            break;
    }
    return Efield;
}