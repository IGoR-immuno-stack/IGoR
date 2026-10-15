/*
 * ScratchDir.h
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <system_error>

/**
 * \brief A scratch directory, removed when the test is done with it.
 *
 * Removal happens in the destructor, so it runs however the test ends, a failed `REQUIRE`
 * included, and never throws: a test that passed does not fail on its cleanup.
 *
 * Declare it before any stream that opens one of its files. Windows refuses to delete a file
 * that is still open, and locals are destroyed in reverse order, so the stream then closes
 * before the directory goes.
 */
class ScratchDir
{
public:
    ScratchDir()
    {
        std::random_device device;
        path_ = std::filesystem::temp_directory_path() / ("igor_test_" + std::to_string(device()));
        std::filesystem::create_directories(path_);
    }
    ~ScratchDir()
    {
        std::error_code ignored;
        std::filesystem::remove_all(path_, ignored);
    }
    ScratchDir(const ScratchDir &) = delete;
    ScratchDir &operator=(const ScratchDir &) = delete;

    /// The path of `name` in the directory, written with `content` first when there is any.
    std::string file(const std::string &name, const std::string &content = "") const
    {
        const std::string full = (path_ / name).string();
        if (!content.empty()) {
            std::ofstream(full) << content;
        }
        return full;
    }

private:
    std::filesystem::path path_;
};
