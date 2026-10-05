/*
 * Platform.h
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

/**
 * \file Platform.h
 * \brief The two operating-system calls IGoR needs under one name on every platform.
 *
 * Out of Legacy/Utils.h since step 1c: the system headers are included here, at global scope,
 * and the shims live in igor::core instead of the global namespace.
 */

#include <cstdint>

#if defined(_WIN32)

#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif

#  include <process.h>
#  include <winsock2.h>
#  include <windows.h>
#  include <ws2tcpip.h>

#else

#  include <unistd.h>

#endif

namespace igor::core {

#if defined(_WIN32)

inline int portable_getpid()
{
    return _getpid();
}

inline uint32_t portable_gethostid()
{
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);

    char hostname[256];
    gethostname(hostname, sizeof(hostname));

    struct addrinfo hints{};
    hints.ai_family = AF_INET;

    struct addrinfo *info;
    if (getaddrinfo(hostname, nullptr, &hints, &info) != 0)
        return 0;

    uint32_t res = ((struct sockaddr_in *)info->ai_addr)->sin_addr.S_un.S_addr;

    freeaddrinfo(info);
    WSACleanup();
    return res;
}

#else

inline int portable_getpid()
{
    return getpid();
}

inline uint32_t portable_gethostid()
{
    return gethostid();
}

#endif

} // namespace igor::core
