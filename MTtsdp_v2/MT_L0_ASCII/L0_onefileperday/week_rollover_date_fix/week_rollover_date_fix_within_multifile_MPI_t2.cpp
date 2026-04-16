#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <regex>
#include <dirent.h> // For directory traversal
#include <sys/stat.h> // For directory existence check
#include <sys/types.h> // For mode_t
#include <libgen.h>
#include <cstring>


// Function to check if a directory exists
bool directory_exists(const std::string& dir) {
    struct stat info;
    if(stat(dir.c_str(), &info) != 0)
        return false;
    return (info.st_mode & S_IFDIR) != 0;
}

// Function to list all .TXT files in a directory
std::vector<std::string> list_txt_files(const std::string& folder) {
    std::vector<std::string> files;
    DIR* dirp = opendir(folder.c_str());
    if (dirp == nullptr) {
        std::cerr << "Failed to open directory: " << folder << std::endl;
        return files;
    }
    struct dirent* dp;
    while ((dp = readdir(dirp)) != nullptr) {
        std::string fname = dp->d_name;
        if (fname.length() >= 4 && fname.substr(fname.length() - 4) == ".TXT") {
            files.push_back(folder + "/" + fname);
        }
    }
    closedir(dirp);
    return files;
}


void ensure_directory_exists(const std::string& dir) {
    struct stat info;
    if (stat(dir.c_str(), &info) != 0) {
        // Directory does not exist, create it
        mkdir(dir.c_str(), 0755);
    } else if (!(info.st_mode & S_IFDIR)) {
        // Path exists but is not a directory
        std::cerr << "Path exists but is not a directory: " << dir << std::endl;
    }
}

void process_file(const std::string& filename) {
    // Extract directory part of filename
    char* filename_cstr = strdup(filename.c_str());
    char* dir_cstr = dirname(filename_cstr);
    std::string dir_path = dir_cstr;
    free(filename_cstr);

    // Ensure the directory exists for output
    ensure_directory_exists(dir_path);

    // Extract the base filename (e.g., 201906160000.TXT)
    filename_cstr = strdup(filename.c_str());
    char* base_cstr = basename(filename_cstr);
    std::string base_name = base_cstr;
    free(filename_cstr);

    // Construct output filename
    std::string output_filename = dir_path + "/processed_" + base_name;

    // Open input file
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Failed to open input file: " << filename << std::endl;
        return;
    }

    // Open output file
    std::ofstream outfile(output_filename);
    if (!outfile) {
        std::cerr << "Failed to open output file: " << output_filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(infile, line)) {
        // Your existing regex extraction and date adjustment code
        std::regex date_regex(R"(^\s*([0-9]+)\s+([0-9]+)\s+([0-9]+))");
        std::smatch match;
        if (std::regex_search(line, match, date_regex)) {
            std::string year = match[1];
            std::string month = match[2];
            std::string day = match[3];

            std::string date_str = year + "-" + month + "-" + day;
            std::string cmd = "date -d \"" + date_str + " +7168 days\" +\"%Y %m %d\"";
            FILE* pipe = popen(cmd.c_str(), "r");
            if (pipe == nullptr) {
                std::cerr << "Failed to run date command" << std::endl;
                continue;
            }
            char buffer[128];
            std::string new_date_str;
            if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
                new_date_str = std::string(buffer);
            }
            pclose(pipe);

            std::istringstream iss(new_date_str);
            std::string new_year, new_month, new_day;
            iss >> new_year >> new_month >> new_day;

            std::string replaced_line = std::regex_replace(line, date_regex,
                new_year + " " + new_month + " " + new_day);

            outfile << replaced_line << "\n";
        } else {
            outfile << line << "\n";
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check if user provided folder names
    std::vector<std::string> folders;
    if (argc > 1) {
	    for (int i = 1; i < argc; ++i) {
            folders.push_back(argv[i]);
        }
    } else {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " folder1 folder2 ...\n";
        }
        MPI_Finalize();
        return 1;
    }

   // Declare all_files here so all processes are aware of it
    std::vector<std::string> all_files;

    if (rank == 0) {
        // Ensure folders exist
        for (const auto& folder : folders) {
            if (!directory_exists(folder)) {
                std::cerr << "Folder does not exist: " << folder << std::endl;
                continue;
            }
            // List all .TXT files in each folder
            auto files_in_folder = list_txt_files(folder);
            all_files.insert(all_files.end(), files_in_folder.begin(), files_in_folder.end());
        }
    }

    // Broadcast total number of files
    int total_files = 0;
    if (rank == 0) {
        total_files = all_files.size();
    }
    MPI_Bcast(&total_files, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Prepare buffer for filenames
    std::vector<int> name_lengths(total_files);
    std::vector<char> all_names_buffer;

    if (rank == 0) {
        for (int i = 0; i < total_files; ++i) {
            name_lengths[i] = all_files[i].size();
        }
        // Concatenate all filenames into one buffer
        int total_length = 0;
        for (const auto& name : all_files) {
            total_length += name.size();
        }
        all_names_buffer.resize(total_length);
        int pos = 0;
        for (const auto& name : all_files) {
            std::copy(name.begin(), name.end(), all_names_buffer.begin() + pos);
            pos += name.size();
        }
    }

    // Broadcast lengths
    MPI_Bcast(name_lengths.data(), total_files, MPI_INT, 0, MPI_COMM_WORLD);
    // Broadcast concatenated buffer size
    int total_buffer_size = 0;
    if (rank != 0) {
        total_buffer_size = 0;
    } else {
        total_buffer_size = all_names_buffer.size();
    }
    MPI_Bcast(&total_buffer_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize buffer on other ranks
    if (rank != 0) {
        all_names_buffer.resize(total_buffer_size);
    }

    // Broadcast the actual filenames buffer
    MPI_Bcast(all_names_buffer.data(), total_buffer_size, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Reconstruct the list of files on all processes
    std::vector<std::string> files;
    size_t offset = 0;
    for (int len : name_lengths) {
        files.push_back(std::string(all_names_buffer.data() + offset, len));
        offset += len;
    }

    // Distribute files among processes
    for (int i = 0; i < total_files; ++i) {
        if (i % size == rank) {
            // Process assigned file
            process_file(files[i]);
            std::cout << "Process " << rank << " processed " << files[i] << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}

