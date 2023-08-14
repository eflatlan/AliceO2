#include "H5Cpp.h"
#include <iostream>
#include <vector>
#include "TH2F.h"
using vecArray2 = std::vector<std::array<double,2>>;
using namespace H5;

class ParticleUtils {
public:




    struct Candidate2 {
        double x, y = 0.;
        int candStatus = 0;
    };

    struct ParticleInfo {

	vecArray2 pionCandidates, kaonCandidates, protonCandidates;
        std::vector<Candidate2> candsCombined;
        std::array<int, 4> arrayInfo;
        double momentum;
        double mass;
        double energy;
        double refractiveIndex;
        double ckov;
        double xRad;
        double yRad;
        double xPC;
        double yPC;
        double thetaP;
        double phiP;
    };

    static void saveParticleInfoToHDF5(std::vector<ParticleInfo>& particleVector) {
        H5File file("ParticleInfo.h5", H5F_ACC_TRUNC);

        H5::CompType mtype(sizeof(Candidate2));
        mtype.insertMember("x", HOFFSET(Candidate2, x), H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("y", HOFFSET(Candidate2, y), H5::PredType::NATIVE_DOUBLE);
        mtype.insertMember("candStatus", HOFFSET(Candidate2, candStatus),  H5::PredType::NATIVE_INT);

				

        for (size_t i = 0; i < particleVector.size(); ++i) {
            auto& particle = particleVector[i];


            std::string groupName = "Particle" + std::to_string(i);

	    // std::cout << groupName << std::endl;

            Group particleGroup = file.createGroup(groupName);

						std::vector<double> x_values;
						std::vector<double> y_values;
						std::vector<int> candStatus_values;
						Printf("ParticleUtils ::: size of candValus = %zu", particle.candsCombined.size());
						for (const auto& cand : particle.candsCombined) {
								x_values.push_back(cand.x);
								y_values.push_back(cand.y);
								candStatus_values.push_back(cand.candStatus);
								//Printf("db x_ %.2f y %.2f cStat %d", cand.x, cand.y, cand.candStatus);
						}
						
						hsize_t dims[1] = { x_values.size() };
						DataSpace dataspace(1, dims);
						DataSet dataset = particleGroup.createDataSet("x_values", PredType::NATIVE_DOUBLE, dataspace);
						dataset.write(x_values.data(), PredType::NATIVE_DOUBLE);

						// Write y_values
						dataspace = DataSpace(1, dims);
						dataset = particleGroup.createDataSet("y_values", PredType::NATIVE_DOUBLE, dataspace);
						dataset.write(&y_values[0], PredType::NATIVE_DOUBLE);

						// Write candStatus_values
						dataspace = DataSpace(1, dims);
						dataset = particleGroup.createDataSet("candStatus_values", PredType::NATIVE_INT, dataspace);
						dataset.write(candStatus_values.data(), PredType::NATIVE_INT);
						
						


            DataSpace attr_dataspace = DataSpace(H5S_SCALAR);

            // Store scalar values
            Attribute attribute = particleGroup.createAttribute("Momentum", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.momentum);

            attribute = particleGroup.createAttribute("Mass", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.mass);

            attribute = particleGroup.createAttribute("Energy", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.energy);

            attribute = particleGroup.createAttribute("RefractiveIndex", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.refractiveIndex);

            attribute = particleGroup.createAttribute("Ckov", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.ckov);

            attribute = particleGroup.createAttribute("xRad", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.xRad);

            attribute = particleGroup.createAttribute("yRad", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.yRad);

            attribute = particleGroup.createAttribute("xPC", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.xPC);

            attribute = particleGroup.createAttribute("yPC", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.yPC);

            attribute = particleGroup.createAttribute("ThetaP", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.thetaP);

            attribute = particleGroup.createAttribute("PhiP", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.phiP);


	    for(int i = 0; i < 4; i++) {
            	std::string str = "ArrayInfo" + std::to_string(i);
            	attribute = particleGroup.createAttribute(str, PredType::NATIVE_DOUBLE, attr_dataspace);
            	attribute.write(PredType::NATIVE_INT, &particle.arrayInfo[i]);
	    }
			/*
	    attribute = particleGroup.createAttribute("C2x", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.candsCombined.x);

	    attribute = particleGroup.createAttribute("C2y", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_DOUBLE, &particle.candsCombined.y);


	    attribute = particleGroup.createAttribute("candStatus", PredType::NATIVE_DOUBLE, attr_dataspace);
            attribute.write(PredType::NATIVE_INT, &particle.candsCombined.candStatus);
	   */ 

	    // this is not working : 
            // Write the candsCombined

            hsize_t vecDims[1] = { particle.candsCombined.size() };
            H5::DataSpace dataspaceVec2(1, vecDims);

            H5::DataSet datasetVec2 = particleGroup.createDataSet("candsCombined", mtype, dataspaceVec2);
            datasetVec2.write(&particle.candsCombined[0], mtype); /**/
            //datasetVec2.write(particle.candsCombined.data(), mtype); /**/


						std::vector<double> X, Y;
						/*for(const auto& c : candsCombined)
						{
							X.emplace_back(c.x);
							Y.emplace_back(c.y);
						}*/


		 				hsize_t vecDimsX[1] = { X.size() };
            H5::DataSpace dataspaceVecX(1, vecDimsX);

            H5::DataSet datasetVecX = particleGroup.createDataSet("candsCombinedX", mtype, dataspaceVecX);
            //datasetVec.write(&particle.candsCombined[0], mtype); /**/
            datasetVecX.write(particle.candsCombined.data(), mtype); /**/





		 				hsize_t vecDimsVec3[1] = { particle.candsCombined.size() };
            H5::DataSpace dataspaceVec3(1, vecDimsVec3);

            H5::DataSet datasetVec3 = particleGroup.createDataSet("candsCombined3", mtype, dataspaceVec3);
            //datasetVec.write(&particle.candsCombined[0], mtype); /**/
            datasetVec3.write(particle.candsCombined.data(), mtype); /**/


		 				hsize_t vecDimsVec4[1] = { particle.candsCombined.size() };
            H5::DataSpace dataspaceVec4(1, vecDimsVec3);

            H5::DataSet datasetVec4 = particleGroup.createDataSet("candsCombined4", mtype, dataspaceVec4);
            datasetVec4.write(&particle.candsCombined[0], mtype); /**/
            //datasetVec3.write(particle.candsCombined.data(), mtype); /**/

        }
    }

    static std::vector<ParticleInfo> loadParticleInfoFromHDF5(const std::string& filename) {
        std::vector<ParticleInfo> particleVector;

        H5File file(filename, H5F_ACC_RDONLY);

        for (hsize_t i = 0; i < file.getNumObjs(); ++i) {
            std::string groupName = file.getObjnameByIdx(i);
            Group particleGroup = file.openGroup(groupName);
						
            // Read simple attributes
            // (The rest of the code is the same as before without the removed parts)

            // Read candsCombined
            H5::CompType mtype(sizeof(Candidate2));
            mtype.insertMember("x", HOFFSET(Candidate2, x), H5::PredType::NATIVE_DOUBLE);
            mtype.insertMember("y", HOFFSET(Candidate2, y), H5::PredType::NATIVE_DOUBLE);
            mtype.insertMember("candStatus", HOFFSET(Candidate2, candStatus), H5::PredType::NATIVE_INT);

            ParticleInfo particleInfo;

            hsize_t dims[1] = { particleInfo.candsCombined.size() };
            H5::DataSpace dataspace(1, dims);
            DataSet dataset = particleGroup.openDataSet("candsCombined");
            dataspace = dataset.getSpace();

            hsize_t dims_out[1];
            dataspace.getSimpleExtentDims(dims_out, NULL);

            std::vector<Candidate2> candsCombined(dims_out[0]);
            dataset.read(candsCombined.data(), mtype);

						for(const auto& c : candsCombined) {Printf("X %.2f Y %.2f cs %d", c.x, c.y, c.candStatus);}

            // (The rest of the code is the same as before without the removed parts)
            //*/ 
            //particleVector.push_back(particleInfo);
        }

        return particleVector;
    }
};
