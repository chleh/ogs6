#include <cstdio>

#include "TESProcess-notpl.h"
#include "TESFEM-notpl.h"


namespace ProcessLib
{

namespace TES
{

void printGlobalMatrix(const Eigen::SparseMatrix<double> &mat)
{
    const unsigned R = mat.rows();
    const unsigned C = mat.cols();
    assert(R==C);

    const unsigned N = R/NODAL_DOF;

    for (unsigned nr=0; nr<N; ++nr)
    {
        for (unsigned vr=0; vr<NODAL_DOF; ++vr)
        {
            const unsigned row = vr*N+nr;

            for (unsigned nc=0; nc<N; ++nc)
            {
                for (unsigned vc=0; vc<NODAL_DOF; ++vc)
                {
                    const unsigned col = vc*N+nc;
                    if (col!=0) std::printf("  ");

                    std::printf("%14g", const_cast<Eigen::SparseMatrix<double>& >(mat).coeffRef(row, col));
                    // std::cout << r << ", " << c << ": " << _A->getRawMatrix().coeffRef(r,c) << std::endl;
                }
            }

            std::printf("\n");
        }
    }


}

} // namespace TES

} // namespace ProcessLib
