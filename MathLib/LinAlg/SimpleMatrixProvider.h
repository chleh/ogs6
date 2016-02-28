// TODO
#pragma once

#include<map>
#include<memory>

#include "MatrixProviderUser.h"

namespace MathLib
{

template<typename Matrix, typename Vector>
class SimpleMatrixProvider final
        : public MatrixProvider<Matrix, Vector>
{
public:
    Matrix& getMatrix(const std::size_t size, std::size_t& id) override;
    void releaseMatrix(std::size_t const id) override;

    Vector& getVector(const std::size_t size, std::size_t& id) override;
    void releaseVector(std::size_t const id);

    void setMatrixSpecificationsProvider(MatrixSpecificationsProvider const& spec_prvd) override;

private:
    std::size_t _next_id = 1;

    std::map<std::size_t, std::unique_ptr<Matrix> > _unused_matrices;
    std::map<std::size_t, std::unique_ptr<Matrix> > _used_matrices;

    std::map<std::size_t, std::unique_ptr<Matrix> > _unused_vectors;
    std::map<std::size_t, std::unique_ptr<Matrix> > _used_vectors;

    MatrixSpecifications _mat_spec = { 0, 0 }; // TODO improve
}



} // namespace MathLib

#include "SimpleMatrixProvider-impl.h"
