// dealii includes
#include <deal.II/numerics/data_out.h>


// Lethe includes
#include <core/manifolds.h>
#include <core/pvd_handler.h>
#include <core/solutions_output.h>


// Std
#include <fstream>
#include <iostream>

template <int dim, int spacedim>
void
write_vtu_and_pvd(PVDHandler                            &pvd_handler,
                  const DataOutInterface<dim, spacedim> &data_out,
                  const std::string                      folder,
                  const std::string                      file_prefix,
                  const double                           time,
                  const unsigned int                     iter,
                  const unsigned int                     group_files,
                  const MPI_Comm                        &mpi_communicator,
                  const unsigned int                     digits)
{
  const int my_id = Utilities::MPI::this_mpi_process(mpi_communicator);

  // Write master files (.pvtu,.pvd,.visit) on the master process
  if (my_id == 0)
    {
      std::vector<std::string> filenames;
      const unsigned int       n_processes =
        Utilities::MPI::n_mpi_processes(mpi_communicator);
      const unsigned int n_files =
        (group_files == 0) ? n_processes : std::min(group_files, n_processes);

      for (unsigned int i = 0; i < n_files; ++i)
        filenames.push_back(file_prefix + "." +
                            Utilities::int_to_string(iter, digits) + "." +
                            Utilities::int_to_string(i, digits) + ".vtu");

      std::string pvtu_filename =
        (file_prefix + "." + Utilities::int_to_string(iter, digits) + ".pvtu");

      std::string   pvtu_filename_with_folder = folder + pvtu_filename;
      std::ofstream master_output(pvtu_filename_with_folder.c_str());

      data_out.write_pvtu_record(master_output, filenames);

      std::string pvdPrefix = (folder + file_prefix + ".pvd");
      pvd_handler.append(time, pvtu_filename);
      std::ofstream pvd_output(pvdPrefix.c_str());
      DataOutBase::write_pvd_record(pvd_output, pvd_handler.times_and_names);
    }

  const unsigned int my_file_id =
    (group_files == 0 ? my_id : my_id % group_files);
  int color = my_id % group_files;

  {
    MPI_Comm comm;
    MPI_Comm_split(mpi_communicator, color, my_id, &comm);
    const std::string filename =
      (folder + file_prefix + "." + Utilities::int_to_string(iter, digits) +
       "." + Utilities::int_to_string(my_file_id, digits) + ".vtu");
    data_out.write_vtu_in_parallel(filename.c_str(), comm);

    MPI_Comm_free(&comm);
  }
}

template <int dim>
void
write_boundaries_vtu(const DataOutFaces<dim> &data_out_faces,
                     const std::string        folder,
                     const double,
                     const unsigned int iter,
                     const MPI_Comm    &mpi_communicator,
                     const std::string  file_prefix,
                     const unsigned int digits)
{
  const int my_id = Utilities::MPI::this_mpi_process(mpi_communicator);

  int      color = my_id % 1;
  MPI_Comm comm;
  MPI_Comm_split(mpi_communicator, color, my_id, &comm);

  const std::string face_filename =
    (folder + file_prefix + "." + Utilities::int_to_string(iter, digits) +
     ".vtu");
  data_out_faces.write_vtu_in_parallel(face_filename.c_str(), comm);

  MPI_Comm_free(&comm);
}

template void
write_vtu_and_pvd(PVDHandler                   &pvd_handler,
                  const DataOutInterface<1, 2> &data_out,
                  const std::string             folder,
                  const std::string             file_prefix,
                  const double                  time,
                  const unsigned int            iter,
                  const unsigned int            group_files,
                  const MPI_Comm               &mpi_communicator,
                  const unsigned int            digits);

template void
write_vtu_and_pvd(PVDHandler                   &pvd_handler,
                  const DataOutInterface<2, 2> &data_out,
                  const std::string             folder,
                  const std::string             file_prefix,
                  const double                  time,
                  const unsigned int            iter,
                  const unsigned int            group_files,
                  const MPI_Comm               &mpi_communicator,
                  const unsigned int            digits);

template void
write_vtu_and_pvd(PVDHandler                   &pvd_handler,
                  const DataOutInterface<2, 3> &data_out,
                  const std::string             folder,
                  const std::string             file_prefix,
                  const double                  time,
                  const unsigned int            iter,
                  const unsigned int            group_files,
                  const MPI_Comm               &mpi_communicator,
                  const unsigned int            digits);

template void
write_vtu_and_pvd(PVDHandler                   &pvd_handler,
                  const DataOutInterface<3, 3> &data_out,
                  const std::string             folder,
                  const std::string             file_prefix,
                  const double                  time,
                  const unsigned int            iter,
                  const unsigned int            group_files,
                  const MPI_Comm               &mpi_communicator,
                  const unsigned int            digits);

template void
write_vtu_and_pvd(PVDHandler                   &pvd_handler,
                  const DataOutInterface<0, 2> &data_out,
                  const std::string             folder,
                  const std::string             file_prefix,
                  const double                  time,
                  const unsigned int            iter,
                  const unsigned int            group_files,
                  const MPI_Comm               &mpi_communicator,
                  const unsigned int            digits);

template void
write_vtu_and_pvd(PVDHandler                   &pvd_handler,
                  const DataOutInterface<0, 3> &data_out,
                  const std::string             folder,
                  const std::string             file_prefix,
                  const double                  time,
                  const unsigned int            iter,
                  const unsigned int            group_files,
                  const MPI_Comm               &mpi_communicator,
                  const unsigned int            digits);


template void
write_boundaries_vtu(const DataOutFaces<2> &data_out_faces,
                     const std::string      folder,
                     const double           time,
                     const unsigned int     iter,
                     const MPI_Comm        &mpi_communicator,
                     const std::string      file_prefix,
                     const unsigned int     digits);

template void
write_boundaries_vtu(const DataOutFaces<3> &data_out_faces,
                     const std::string      folder,
                     const double           time,
                     const unsigned int     iter,
                     const MPI_Comm        &mpi_communicator,
                     const std::string      file_prefix,
                     const unsigned int     digits);
