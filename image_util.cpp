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

#include <string>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <stdexcept>
#include "image_util.h"

JPGCompressor::JPGCompressor(int quality, bool flipy)
  : compressor(tjInitCompress()), buffer(nullptr), bufsize(0), quality(quality), flipy(flipy)
{}
JPGCompressor::~JPGCompressor() {
  if (buffer) {
    tjFree(buffer);
  }
  tjDestroy(compressor);
}
const std::pair<unsigned char*, unsigned long> JPGCompressor::compress(uint32_t *pixels,
    int width, int height)
{
  const int flags = flipy ? TJFLAG_BOTTOMUP : 0;
  const int rc = tjCompress2(compressor, reinterpret_cast<unsigned char*>(pixels),
      width, width * 4, height, TJPF_RGBA, &buffer, &bufsize, TJSAMP_420,
      quality, flags);
  if (rc != 0) {
    const std::string tj_err = tjGetErrorStr();
    throw std::runtime_error("Failed to compress JPG! Error: " + tj_err);
  }
  return std::make_pair(buffer, bufsize);
}

JPGDecompressor::JPGDecompressor() : decompressor(tjInitDecompress()) {}
JPGDecompressor::~JPGDecompressor() {
  tjDestroy(decompressor);
}
void JPGDecompressor::decompress(unsigned char *jpg, const unsigned long jpeg_size,
    const int width, const int height, std::vector<uint32_t> &img)
{
  if (img.size() != width * height * 4) {
    img.resize(width * height * 4, 0);
  }
  const int rc = tjDecompress2(decompressor, jpg, jpeg_size,
      reinterpret_cast<unsigned char*>(img.data()),
      width, width * 4, height, TJPF_RGBA, TJFLAG_BOTTOMUP);
  if (rc != 0) {
    const std::string tj_err = tjGetErrorStr();
    throw std::runtime_error("Failed to decompress JPG! Error: " + tj_err);
  }
}

void save_jpeg_file(const std::string &fname, uint32_t *pixels,
    int width, int height)
{
  JPGCompressor compressor(100);
  auto jpg = compressor.compress(pixels, width, height);
  std::ofstream fout(fname.c_str(), std::ios::binary);
  fout.write(reinterpret_cast<const char*>(jpg.first), jpg.second);
}


