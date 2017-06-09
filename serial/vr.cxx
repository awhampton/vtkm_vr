#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>

#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/DeviceAdapterAlgorithm.h>

#include <vtkm/Math.h>
#include <vtkm/VectorAnalysis.h>

#include <iostream>

#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkDataSetWriter.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSphereSource.h>

#include <vtkCamera.h>
#include <vtkDataSetMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>




// ****************************************************************************
//                           GLOBAL VARIABLES
// ****************************************************************************

// these are set as command line arguments
const vtkm::Int32 NUM_SAMPLES = 1000;                // e.g., 100
const vtkm::Int32 WIDTH = 1024;                      // e.g., 1024
const vtkm::Int32 HEIGHT = 1024;                     // e.g., 1024
const std::string filename = "../../astro512.vtk";      // e.g., "../astro64.vtk"
const std::string out_filename = "nova.png";         // e.g., "nova.png";

// default samples affects the color intensity
// smaller default samples --> pixel saturated more quickly
const vtkm::Int32    DEFAULT_SAMPLES = 250;  // 250 for ASTRO is good
const vtkm::Float64  SAMPLE_RATIO = ((vtkm::Float64) DEFAULT_SAMPLES) / NUM_SAMPLES;



// ****************************************************************************
//                         CAMERA & TRANSFER FUNCTION
// ****************************************************************************

class Camera
{
public:
    vtkm::Float64      near, far;
    vtkm::Float64      angle;
    vtkm::Vec<vtkm::Float64, 3>  position;
    vtkm::Vec<vtkm::Float64, 3>  focus;
    vtkm::Vec<vtkm::Float64, 3>  up;

    VTKM_CONT
    Camera(){
      focus[0] = 0;
      focus[1] = 0;
      focus[2] = 0;
      up[0] = 0;
      up[1] = 1;
      up[2] = 0;
      angle = 30;
      near = 7.5e+7;
      far = 1.4e+8;
      position[0] = -8.25e+7;
      position[1] = -3.45e+7;
      position[2] = 3.35e+7;
    }
};

class TransferFunction
{
public:

    vtkm::Float64      min;
    vtkm::Float64      max;
    vtkm::Int32        numBins;
    vtkm::UInt8        colors[3*256];    // size is 3*numBins
    vtkm::Float64      opacities[256]; // size is numBins
 
    VTKM_CONT
    TransferFunction(){

      vtkm::Int32  i;

      min = 10;
      max = 15;
      numBins = 256;
      vtkm::UInt8 charOpacity[256] = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 17, 17, 17, 17, 17, 17, 16, 16, 15, 14, 13, 12, 11, 9, 8, 7, 6, 5, 5, 4, 3, 3, 3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 16, 18, 20, 22, 24, 27, 29, 32, 35, 38, 41, 44, 47, 50, 52, 55, 58, 60, 62, 64, 66, 67, 68, 69, 70, 70, 70, 69, 68, 67, 66, 64, 62, 60, 58, 55, 52, 50, 47, 44, 41, 38, 35, 32, 29, 27, 24, 22, 20, 20, 23, 28, 33, 38, 45, 51, 59, 67, 76, 85, 95, 105, 116, 127, 138, 149, 160, 170, 180, 189, 198, 205, 212, 217, 221, 223, 224, 224, 222, 219, 214, 208, 201, 193, 184, 174, 164, 153, 142, 131, 120, 109, 99, 89, 79, 70, 62, 54, 47, 40, 35, 30, 25, 21, 17, 14, 12, 10, 8, 6, 5, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        };

      for (i = 0 ; i < 256 ; i++)
        opacities[i] = charOpacity[i]/255.0;
      const vtkm::Int32 numControlPoints = 8;
      vtkm::UInt8 controlPointColors[numControlPoints*3] = {
           71, 71, 219, 0, 0, 91, 0, 255, 255, 0, 127, 0,
           255, 255, 0, 255, 96, 0, 107, 0, 0, 224, 76, 76
       };
      vtkm::Float64 controlPointPositions[numControlPoints] = { 0, 0.143, 0.285, 0.429, 0.571, 0.714, 0.857, 1.0 };
      for (i = 0 ; i < numControlPoints-1 ; i++)
      {
        vtkm::Int32 start = controlPointPositions[i]*numBins;
        vtkm::Int32 end   = controlPointPositions[i+1]*numBins+1;
        if (end >= numBins)
            end = numBins-1;
        for (vtkm::Int32 j = start ; j <= end ; j++)
        {
            vtkm::Float64 proportion = (j/(numBins-1.0)-controlPointPositions[i])/(controlPointPositions[i+1]-controlPointPositions[i]);
            if (proportion < 0 || proportion > 1.)
                continue;
            for (vtkm::Int32 k = 0 ; k < 3 ; k++)
                colors[3*j+k] = proportion*(controlPointColors[3*(i+1)+k]-controlPointColors[3*i+k])
                                 + controlPointColors[3*i+k];
        }
    }
  }

    VTKM_EXEC
    vtkm::Int32 GetBin(vtkm::Float64 value) const
    {
        vtkm::Float64 step = (max - min) / numBins;
        return (vtkm::Int32) ((value - min) / step);
    }

    VTKM_EXEC
    void ApplyTransferFunction(vtkm::Float64 value, vtkm::Vec<vtkm::UInt8, 3> &RGB, vtkm::Float64 &opacity) const
    {
        if(value < min || value > max)
        {
            RGB[0] = 0;
            RGB[1] = 0;
            RGB[2] = 0;
            opacity = 0;
        }
        else
        {
            vtkm::Int32 bin = GetBin(value);
            #ifdef DEBUG
                std::cout << "tf bin: " << bin << std::endl;
            #endif
            RGB[0] = colors[3*bin+0];
            RGB[1] = colors[3*bin+1];
            RGB[2] = colors[3*bin+2];
            opacity = opacities[bin];
        }
    }
};



// ****************************************************************************
//                         VOLUME RENDERER WORKLET
// ****************************************************************************

namespace vtkm{
namespace worklet{

class VolumeRenderer : public vtkm::worklet::WorkletMapField
{
public:
    typedef void ControlSignature(FieldIn<> v1, WholeArrayIn<> X, WholeArrayIn<> Y, WholeArrayIn<> Z, WholeArrayIn<> F, FieldOut<> result);
    typedef void ExecutionSignature(_1, _2, _3, _4, _5, _6, WorkIndex);
    typedef _1 InputDomain;
    TransferFunction tf;
    Camera cam;
    vtkm::Vec<vtkm::Int32, 3> dims;
    vtkm::Float32 X_min, X_max, Y_min, Y_max, Z_min, Z_max;
    vtkm::Vec<vtkm::Float64, 3> look;
    vtkm::Vec<vtkm::Float64, 3> u;
    vtkm::Vec<vtkm::Float64, 3> v;
    vtkm::Vec<vtkm::Float64, 3> d_x;
    vtkm::Vec<vtkm::Float64, 3> d_y;

    // constructor
    VTKM_CONT
    VolumeRenderer(vtkm::Vec<vtkm::Int32, 3> in_dims, vtkm::Int32 x1, vtkm::Int32 x2, vtkm::Int32 y1, vtkm::Int32 y2, vtkm::Int32 z1, vtkm::Int32 z2){

        // store the dataset dimensions and min/max values
        dims[0] = in_dims[0];
        dims[1] = in_dims[1];
        dims[2] = in_dims[2];

        X_min = x1;
        X_max = x2;

        Y_min = y1;
        Y_max = y2;

        Z_min = z1;
        Z_max = z2;

        // pre-compute some of the vectors for finding each ray
        look = cam.focus - cam.position;
        vtkm::Normalize(look);

        u = vtkm::Cross(look, cam.up);
        vtkm::Normalize(u);

        v = vtkm::Cross(look, u);
        vtkm::Normalize(v);

        vtkm::Float64 d_x_coeff = (tan(cam.angle * 3.1415 / (2 * 180.0))) / ((vtkm::Float64) WIDTH);
        d_x = d_x_coeff * u;

        vtkm::Float64 d_y_coeff = (tan(cam.angle * 3.1415 / (2 * 180.0))) / ((vtkm::Float64) HEIGHT);
        d_y = d_y_coeff * v;
    }



    // ****************************************************************************
    //                           INTERPOLATE VALUES
    // ****************************************************************************

    // note: if a == b, return the value at a
    VTKM_EXEC
    vtkm::Float64 InterpolateValues(vtkm::Float64 A, vtkm::Float64 B, vtkm::Float64 a, vtkm::Float64 b, vtkm::Float64 x) const
    {
        if( b - a < 0.00001 )
        {
            return A;
        }
        else
        {
            vtkm::Float64 t = (x - a) / (b - a);
            return A + t * (B - A);
        }
    }



    // ****************************************************************************
    //                           DATASET FUNCTIONS
    // ****************************************************************************

    // ****************************************************************************
    //  Function: GetNumberOfPoints
    //
    //  Arguments:
    //     dims: an array of size 3 with the number of points in X, Y, and Z.
    //           2D data sets would have Z=1
    //
    //  Returns:  the number of points in a rectilinear mesh
    //
    // ****************************************************************************

    VTKM_EXEC
    vtkm::Int32 GetNumberOfPoints(const vtkm::Vec<vtkm::Int32, 3> dims) const
    {
        // 3D
        return dims[0]*dims[1]*dims[2];
    }

    // ****************************************************************************
    //  Function: GetNumberOfCells
    //
    //  Arguments:
    //
    //      dims: an array of size 3 with the number of points in X, Y, and Z.
    //            2D data sets would have Z=1
    //
    //  Returns:  the number of cells in a rectilinear mesh
    //
    // ****************************************************************************

    VTKM_EXEC
    vtkm::Int32 GetNumberOfCells(const vtkm::Vec<vtkm::Int32, 3> dims) const
    {
        // 3D
        return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    }

    // ****************************************************************************
    //  Function: GetPointIndex
    //
    //  Arguments:
    //      idx:  the logical index of a point.
    //              0 <= idx[0] < dims[0]
    //              1 <= idx[1] < dims[1]
    //              2 <= idx[2] < dims[2] (or always 0 if 2D)
    //      dims: an array of size 3 with the number of points in X, Y, and Z.
    //            2D data sets would have Z=1
    //
    //  Returns:  the point index
    //
    // ****************************************************************************

    VTKM_EXEC
    vtkm::Int32 GetPointIndex(const vtkm::Int32 *idx, const vtkm::Vec<vtkm::Int32, 3> dims) const
    {
        // 3D
        return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    }

    // ****************************************************************************
    //  Function: GetCellIndex
    //
    //  Arguments:
    //      idx:  the logical index of a cell.
    //              0 <= idx[0] < dims[0]-1
    //              1 <= idx[1] < dims[1]-1
    //              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
    //      dims: an array of size 3 with the number of points in X, Y, and Z.
    //            2D data sets would have Z=1
    //
    //  Returns:  the cell index
    //
    // ****************************************************************************

    VTKM_EXEC
    vtkm::Int32 GetCellIndex(const vtkm::Int32 *idx, const vtkm::Vec<vtkm::Int32, 3> dims) const
    {
        // 3D
        return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    }

    // ****************************************************************************
    //  Function: GetLogicalPointIndex
    //
    //  Arguments:
    //      idx (output):  the logical index of the point.
    //              0 <= idx[0] < dims[0]
    //              1 <= idx[1] < dims[1]
    //              2 <= idx[2] < dims[2] (or always 0 if 2D)
    //      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
    //      dims: an array of size 3 with the number of points in X, Y, and Z.
    //            2D data sets would have Z=1
    //
    //  Returns:  None (argument idx is output)
    //
    // ****************************************************************************

    VTKM_EXEC
    void GetLogicalPointIndex(vtkm::Int32 *idx, vtkm::Int32 pointId, const vtkm::Vec<vtkm::Int32,3 > dims) const
    {
        // 3D
        idx[0] = pointId%dims[0];
        idx[1] = (pointId/dims[0])%dims[1];
        idx[2] = pointId/(dims[0]*dims[1]);

    }

    // ****************************************************************************
    //  Function: GetLogicalCellIndex
    //
    //  Arguments:
    //      idx (output):  the logical index of the cell index.
    //              0 <= idx[0] < dims[0]-1
    //              1 <= idx[1] < dims[1]-1
    //              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
    //      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
    //      dims: an array of size 3 with the number of points in X, Y, and Z.
    //            2D data sets would have Z=1
    //
    //  Returns:  None (argument idx is output)
    //
    // ****************************************************************************

    VTKM_EXEC
    void GetLogicalCellIndex(vtkm::Int32 *idx, vtkm::Int32 cellId, const vtkm::Vec<vtkm::Int32, 3> dims) const
    {
        // 3D
        idx[0] = cellId%(dims[0]-1);
        idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
        idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    }

    // each dimension should be at least 2
    //   note: assumes that pt is actually contained in some cell
    //         this check should be performed separately
    template<typename InFieldType>
    VTKM_EXEC
    vtkm::Float32
    EvaluateFieldAtLocation(const vtkm::Vec<vtkm::Float64, 3> pt, const vtkm::Vec<vtkm::Int32, 3> dims,
                            const InFieldType X, const InFieldType Y,  const InFieldType Z, const InFieldType F) const
    {
        // find the cell that pt falls in
        vtkm::Int32 x_idx = 0;
        vtkm::Int32 y_idx = 0;
        vtkm::Int32 z_idx = 0;

        while( x_idx < dims[0] - 1 && X.Get(x_idx + 1) <= pt[0] )
        {
            x_idx++;
        }
        while( y_idx < dims[1] - 1 && Y.Get(y_idx + 1) <= pt[1] )
        {
            y_idx++;
        }
        while( z_idx < dims[2] - 1 && Z.Get(z_idx + 1) <= pt[2] )
        {
            z_idx++;
        }

        // interpolate the field!
        vtkm::Vec<vtkm::Float32, 8> fbox;
        vtkm::Int32 tmp_idx[3];
        for(vtkm::Int32 xo = 0; xo <= 1; xo++){
          for(vtkm::Int32 yo = 0; yo <= 1; yo++){
            for(vtkm::Int32 zo = 0; zo <= 1; zo++){
              tmp_idx[0] = x_idx + xo;
              tmp_idx[1] = y_idx + yo;
              tmp_idx[2] = z_idx + zo;
              fbox[4 * xo + 2 * yo + zo] = F.Get(GetPointIndex(tmp_idx, dims));
            }
          }
        }

        vtkm::Float64 edge1 = InterpolateValues(fbox[0], fbox[4], X.Get(x_idx), X.Get(x_idx + 1), pt[0]);
        vtkm::Float64 edge2 = InterpolateValues(fbox[1], fbox[5], X.Get(x_idx), X.Get(x_idx + 1), pt[0]);
        vtkm::Float64 edge3 = InterpolateValues(fbox[2], fbox[6], X.Get(x_idx), X.Get(x_idx + 1), pt[0]);
        vtkm::Float64 edge4 = InterpolateValues(fbox[3], fbox[7], X.Get(x_idx), X.Get(x_idx + 1), pt[0]);

        vtkm::Float64 face1 = InterpolateValues(edge1, edge3, Y.Get(y_idx), Y.Get(y_idx + 1), pt[1]);
        vtkm::Float64 face2 = InterpolateValues(edge2, edge4, Y.Get(y_idx), Y.Get(y_idx + 1), pt[1]);

        return InterpolateValues(face1, face2, Z.Get(z_idx), Z.Get(z_idx + 1), pt[2]);
    }



    // ****************************************************************************
    //                           RENDERING COMPUTATIONS
    // ****************************************************************************

    // calculates the ray through pixel (x,y)
    VTKM_EXEC
    void
    CalculateRay(vtkm::Int32 x, vtkm::Int32 y, vtkm::Int32 W, vtkm::Int32 H, vtkm::Vec<vtkm::Float64, 3> &res) const
    {
        res = look + (2*x + 1 - W) * d_x + (2*y + 1 - H) * d_y;
        vtkm::Normalize(res);
    }

    // takes samples in the direction of r, starting at the camera's position
    // between two spheres of radius camera.near and camera.far
    VTKM_EXEC
    void
    GetSampleLocations(Camera c, vtkm::Vec<vtkm::Float64, 3> r, vtkm::Vec<vtkm::Float64, 3> *res) const
    {
        vtkm::Float64 step = (c.far - c.near) / NUM_SAMPLES;
        res[0] = c.position + c.near * r;
        for (vtkm::Int32 i = 1; i < NUM_SAMPLES; i++)
        {
            res[i] = res[i-1] + step * r;
        }
    }

    // performs the sampling through a pixel
    template<typename InFieldType, typename OutFieldType>
    VTKM_EXEC
    void operator() (const OutFieldType &in_values, const InFieldType &X, const InFieldType &Y, const InFieldType &Z, const InFieldType &F, OutFieldType &result, vtkm::Id workIndex) const
    {

        /* don't loop! instead recover (i,j) from thread id */
        //   warning: VTK writes the first row on the bottom of the image!
        vtkm::Int32 j = HEIGHT - 1 - (workIndex / WIDTH);  // row
        vtkm::Int32 i = workIndex % WIDTH;                 // column

        // debug info
        #ifdef DEBUG
            std::cout << std::endl;
            std::cout << "pixel: " << i << " " << j << std::endl;
        #endif

        // calculate ray
        vtkm::Vec<vtkm::Float64, 3> ray;
        CalculateRay(i, j, WIDTH, HEIGHT, ray);
        #ifdef DEBUG
            std::cout << "ray: " << ray[0] << " " << ray[1] << " " << ray[2] << std::endl;
            std::cout << "=====================================================" << std::endl;
        #endif

        // calculate sample locations along ray
        vtkm::Vec<vtkm::Float64, 3> sample_locations[NUM_SAMPLES];
        GetSampleLocations(cam, ray, sample_locations);

        // intersect volume with ray
        vtkm::Vec<vtkm::UInt8, 3> RGB = vtkm::make_Vec(0, 0, 0);
        vtkm::Vec<vtkm::Int32, 3> white = vtkm::make_Vec(255, 255, 255);

        vtkm::Float64 opacity = 0;
        for(vtkm::Int32 k = 0; k < NUM_SAMPLES; k++)
        {

            // get the sample location
            vtkm::Vec<vtkm::Float64, 3> sample = sample_locations[k];
            #ifdef DEBUG
                std::cout << "sample locations: " << sample[0] << " " << sample[1] << " " << sample[2] << std::endl;
            #endif

            // check if the sample location is within the dataset
            if(sample[0] < X_min || sample[0] > X_max || sample[1] < Y_min || sample[1] > Y_max || sample[2] < Z_min || sample[2] > Z_max)
            {
                #ifdef DEBUG
                    std::cout << "skipping" << std::endl;
                #endif
                continue;
            }

            // get the field value
            vtkm::Float64 field_val = EvaluateFieldAtLocation(sample, dims, X, Y, Z, F);
            #ifdef DEBUG
                std::cout << "field value: " << field_val << std::endl;
            #endif

            // apply the transfer function
            vtkm::Vec<vtkm::UInt8, 3> RGB_new;
            vtkm::Float64 opacity_new;
            tf.ApplyTransferFunction(field_val, RGB_new, opacity_new);

            // apply the opacity correction
            opacity_new = 1 - vtkm::Pow((1 - opacity_new), SAMPLE_RATIO);
            #ifdef DEBUG
                std::cout << "corrected opacity: " << opacity_new << std::endl;
            #endif

            // blend with existing value
            RGB = vtkm::Min(white, (vtkm::Vec<vtkm::Int32, 3>) RGB + (1 - opacity) * opacity_new * RGB_new);
            opacity += (1 - opacity) * opacity_new;

            // break early if pixel is totally saturated
            if(opacity > 0.999){
                break;
            }
        }

        // store in result
        result = RGB;


    }
};

}
}



// ****************************************************************************
//                           VTK IMAGE FUNCTIONS
// ****************************************************************************

void
WriteImage(vtkImageData *img, std::string filename)
{
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(filename.c_str());
    writer->Write();
    writer->Delete();
}


vtkImageData *
NewImage(vtkm::Int32 width, vtkm::Int32 height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return image;
}



// ****************************************************************************
//                           MAIN
// ****************************************************************************

int main(int argc, char **argv)
{
    // set up the image
    vtkImageData *image = NewImage(WIDTH, HEIGHT);
    vtkm::UInt8 *buffer = (vtkm::UInt8 *) image->GetScalarPointer(0,0,0);

    // read the input data
    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName(filename.c_str());
    rdr->Update();

    // read the input data as VTK-m dataset
    // TODO: try to use VTK-m datasets for the whole thing
    //vtkm::io::reader::VTKDataSetReader vtkm_reader(filename.c_str());
    //vtkm::cont::DataSet inData = vtkm_reader.ReadDataSet();
    //inData.PrintSummary(std::cout);

    // read the dataset dimensions
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    vtkm::Int32 dims[3];
    rgrid->GetDimensions(dims);
    vtkm::Vec<vtkm::Int32, 3> input_dims(dims[0], dims[1], dims[2]);

    // debug output: range of scalar field
    #ifdef DEBUG
        std::cout << "scalars range: " << (rgrid->GetPointData()->GetScalars()->GetRange())[0] << " -- " << (rgrid->GetPointData()->GetScalars()->GetRange())[1] << std::endl;
    #endif

    // extract the cells and field
    vtkm::Float32 *X = (vtkm::Float32 *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    vtkm::Float32 *Y = (vtkm::Float32 *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    vtkm::Float32 *Z = (vtkm::Float32 *) rgrid->GetZCoordinates()->GetVoidPointer(0);
    vtkm::Float32 *F = (vtkm::Float32 *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    vtkm::Int32 n_field = rgrid->GetNumberOfPoints();

    // get the max and min values
    vtkm::Int32 X_min = X[0];
    vtkm::Int32 X_max = X[dims[0] - 1];

    vtkm::Int32 Y_min = Y[0];
    vtkm::Int32 Y_max = Y[dims[1] - 1];

    vtkm::Int32 Z_min = Z[0];
    vtkm::Int32 Z_max = Z[dims[2] - 1];

    // make array handles
    vtkm::Vec<vtkm::UInt8, 3> in_buffer[WIDTH * HEIGHT];
    vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::UInt8, 3> > in_buffer_ah = vtkm::cont::make_ArrayHandle(in_buffer, WIDTH * HEIGHT);
    vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::UInt8, 3> > out_buffer_ah;
    vtkm::cont::ArrayHandle< vtkm::Float32 > X_ah = vtkm::cont::make_ArrayHandle(X, input_dims[0]);
    vtkm::cont::ArrayHandle< vtkm::Float32 > Y_ah = vtkm::cont::make_ArrayHandle(Y, input_dims[1]);
    vtkm::cont::ArrayHandle< vtkm::Float32 > Z_ah = vtkm::cont::make_ArrayHandle(Z, input_dims[2]);
    vtkm::cont::ArrayHandle< vtkm::Float32 > F_ah = vtkm::cont::make_ArrayHandle(F, n_field);

    // instantiate and invoke the volume renderer
    vtkm::worklet::VolumeRenderer vr(input_dims, X_min, X_max, Y_min, Y_max, Z_min, Z_max);
    vtkm::worklet::DispatcherMapField<vtkm::worklet::VolumeRenderer> dispatcher(vr);
    dispatcher.Invoke(in_buffer_ah, X_ah, Y_ah, Z_ah, F_ah, out_buffer_ah);

    // dump the output into buffer
    // TODO: try to use buffer right in the algorithm to avoid copy
    typedef vtkm::cont::ArrayHandle< vtkm::Vec<vtkm::UInt8, 3> >::PortalControl portalType;
    portalType portal = out_buffer_ah.GetPortalControl();
    for(vtkm::Id idx = 0; idx < portal.GetNumberOfValues(); idx++){
      vtkm::Vec<vtkm::UInt8, 3> v = portal.Get(idx);
      buffer[3 * idx + 0] = v[0];
      buffer[3 * idx + 1] = v[1];
      buffer[3 * idx + 2] = v[2];
    }

    // write the image
    WriteImage(image, out_filename);

    return 0;
}
