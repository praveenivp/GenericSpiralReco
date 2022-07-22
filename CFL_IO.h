#pragma once 

#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include <fstream>
#include <complex>
#include <utility>


//only required for hoNDarray
#include "hoNDArray.h"
#define USE_HONDARRAY 1


namespace CFL_IO
{

    bool writeHeader(std::string & filename, std::vector<uint16_t>& dims)
    {
        std::ofstream ofile;
        std::string header_fn = filename + std::string(".hdr");
        ofile.open(header_fn, std::ofstream::trunc);
        ofile << "# Dimensions\n";
        struct remove_comma : std::numpunct<char> {
         char do_thousands_sep()   const { return '\0'; }  // separate with spaces
        };

        ofile.imbue(std::locale(ofile.getloc(), new remove_comma));
        for (auto cdim : dims)
            ofile << cdim<< ' ';
        for (int i = 0; i < 10 - dims.size(); i++)
            ofile << 1 << ' ';
        ofile << std::endl;
        ofile.close();
        return true;
    }

    std::vector<uint16_t> readHeader(std::string &filename)
    {
        std::vector<uint16_t> dims(10, 1);
        std::fstream ifile;
        std::string header_fn = filename + std::string(".hdr");
        ifile.open(header_fn);
        std::string temp;
        getline(ifile, temp, '\n');
        if (strcmp(temp.c_str(), "# Dimensions\n"))
        {
            ifile >> dims[0] >> dims[1] >> dims[2] >> dims[3] >> dims[4] >> dims[5] >> dims[6] >> dims[7] >> dims[8] >> dims[9];
            getline(ifile, temp, '\n');
        }
        else
            std::cerr << "error parsing" << std::endl;

        ifile.close();

        return dims;
    }

    template <typename T>
    bool writeCFL(std::string filename, std::vector<T>& data, std::vector<uint16_t>& dims)
    {
        CFL_IO::writeHeader(filename, dims);
        std::ofstream ofile;
        std::string header_fn = filename + std::string(".cfl");
        ofile.open(header_fn, std::ofstream::trunc | std::ofstream::binary);
        ofile.write(reinterpret_cast<char*>(&data[0]),sizeof(T)*data.size());
        ofile.close();
        return true;
    }

    template <typename T>
    std::pair<std::vector<std::complex<float>>, std::vector<uint16_t>> readCFL(std::string filename)
    {
        auto dims = CFL_IO::readHeader(filename);
        uint64_t numElem = 1;
        for (auto cdim : dims)
            numElem *= cdim;
        
        std::vector<T> cfl_buff(numElem,std::complex(0.0f,0.0f));

        std::ifstream ifile;
        std::string cfl_fn = filename + std::string(".cfl");
        ifile.open(cfl_fn,std::ifstream::binary);
        ifile.read(reinterpret_cast<char*>(&cfl_buff[0]),sizeof(T)*numElem);
        ifile.close();

        return {cfl_buff,dims};
    }



    #ifdef USE_HONDARRAY

            template <typename T>
    bool hoNDArray2CFL(std::string filename, Gadgetron::hoNDArray<T> & data)
    {
        auto Ndims =data.get_number_of_dimensions();
        std::vector<u_int16_t> dims(Ndims,0);
        for (auto i=0; i<Ndims;i++)
            dims[i]=data.get_size(i);
        CFL_IO::writeHeader(filename, dims);
        std::ofstream ofile;
        std::string cfl_fn = filename + std::string(".cfl");
        ofile.open(cfl_fn, std::ofstream::trunc | std::ofstream::binary);
        ofile.write(reinterpret_cast<char*>(data.get_data_ptr()),sizeof(T)*data.get_number_of_elements());
        ofile.close();
        return true;
    }

    template <typename T>
    Gadgetron::hoNDArray<T> CFL2hoNDARRAY(std::string filename)
    {
        auto dims = CFL_IO::readHeader(filename);
        uint64_t numElem = 1;
        for (auto cdim : dims)
        {
            numElem *= cdim;
        }
        
        Gadgetron::hoNDArray<T> cfl_buff(dims[0],dims[1],dims[2],dims[3],dims[4],dims[5],dims[6]);

        std::ifstream ifile;
        std::string cfl_fn = filename + std::string(".cfl");
        ifile.open(cfl_fn,std::ifstream::binary);
        ifile.read(reinterpret_cast<char*>(&cfl_buff[0]),sizeof(T)*numElem);
        ifile.close();
        return cfl_buff;
    }

    #endif

};

// int main()
// {
//     std::string fn("/media/sf_VM_Shared/gadgetronstuff/testCFL");
//     std::complex<float> *ptr = NULL;

//     // create test data
//     std::vector<uint16_t> dims{4, 2, 1};
//     uint64_t numElem = 1;
//     for (auto cdim : dims)
//         numElem *= cdim;

//     std::vector<std::complex<float>> testData(numElem,std::complex(0.0f,0.0f));
//     for (uint64_t i = 0; i < numElem; ++i)
//     {
//         testData[i].real(0.1f * (float)i);
//         testData[i].imag(0.2 * (float)i);
//     }

//     CFL_IO::writeCFL<std::complex<float>>(fn, testData, dims);

//     auto test =CFL_IO::readCFL<std::complex<float>>(fn);
//     auto cfl_data = test.first;
//     auto cfl_hdr = test.first;

//     std::cout<<" cfl data"<<std::endl;
//     for (int i=0;i<cfl_data.size();i++)
//         std::cout<<cfl_data[i].real()<<"+i"<<cfl_data[i].imag()<< " , ";
//     std::cout<<std::endl;

//     std::cout<<" cfl header"<<std::endl;
//             for (auto tt : dims)
//             std::cout << tt << " ";
//         std::cout << std::endl;
//     return 0;
// }
