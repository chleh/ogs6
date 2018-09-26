/** This file is taken from Oliver Sander's book draft! It must not be published
 * as is!
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif  // HAVE_CONFIG_H

#include "dune-config.h"

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

using namespace Dune;

int main(int argc, char* argv[]) try
{
    // Set up MPI, if available
    MPIHelper::instance(argc, argv);

    // ////////////////////////////////
    //   Generate the grid
    // ////////////////////////////////

    if (argc != 3)
    {
        std::cerr << "exactly two arguments taken\n";
        return 1;
    }
    std::string const infile = argv[1];
    std::string const outfile_no_ext = argv[2];
    std::cout << "reading mesh " << infile << '\n';

    const int dim = 2;
    typedef UGGrid<dim> GridType;
    std::shared_ptr<GridType> grid(GmshReader<GridType>::read(infile));

    // grid->globalRefine(2);

    typedef GridType::LeafGridView GridView;
    GridView gridView = grid->leafGridView();

    // Output result
    std::cout << "writing mesh " << outfile_no_ext << '\n';
    VTKWriter<GridView> vtkWriter(gridView);
    // vtkWriter.addVertexData(x, "solution");
    vtkWriter.write(outfile_no_ext);

    return 0;
}
// Error handling
catch (std::exception& e)
{
    std::cout << e.what() << std::endl;
    return 1;
}
