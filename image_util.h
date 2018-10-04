// ======================================================================== //
// Copyright 2018 Intel Corporation                                         //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

#pragma once

#include <utility>
#include <vector>
#include <turbojpeg.h>

class JPGCompressor {
  tjhandle compressor;
  unsigned char *buffer;
  unsigned long bufsize;
  int quality;
  bool flipy;

public:
  JPGCompressor(int quality, bool flipy = true);
  ~JPGCompressor();

  /* Compress and RGBA image and return a pointer to the JPG buffer.
   * The pointer will be valid for the image until the next time compress is
   * called, as the buffer will be re-used.
   */
  const std::pair<unsigned char*, unsigned long> compress(uint32_t *pixels,
      int width, int height);
};

class JPGDecompressor {
  tjhandle decompressor;

public:
  JPGDecompressor();
  ~JPGDecompressor();

  /* Decompress an RGBA JPG image into the buffer passed. img may be resized
   * to ensure it has enough room.
   */
  void decompress(unsigned char *jpg, const unsigned long jpeg_size,
    const int width, const int height, std::vector<uint32_t> &img);
};

void save_jpeg_file(const std::string &fname, uint32_t *pixels,
    int width, int height);

