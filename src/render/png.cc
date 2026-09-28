#include "render/png.h"

#include <zlib.h>

#include <cstring>
#include <fstream>
#include <iterator>
#include <stdexcept>

namespace ssp {
namespace {

constexpr uint8_t kSignature[8] = {0x89, 'P', 'N', 'G', '\r', '\n', 0x1a, '\n'};

void put32(std::vector<uint8_t>& out, uint32_t v) {
    for (int s = 24; s >= 0; s -= 8) out.push_back(static_cast<uint8_t>(v >> s));
}

uint32_t get32(const uint8_t* p) {
    return (uint32_t{p[0]} << 24) | (uint32_t{p[1]} << 16) | (uint32_t{p[2]} << 8) | p[3];
}

void put_chunk(std::vector<uint8_t>& out, const char* type, const std::vector<uint8_t>& body) {
    put32(out, static_cast<uint32_t>(body.size()));
    const std::size_t start = out.size();
    out.insert(out.end(), type, type + 4);
    out.insert(out.end(), body.begin(), body.end());
    put32(out, static_cast<uint32_t>(crc32(0, out.data() + start, out.size() - start)));
}

}  // namespace

std::vector<uint8_t> encode_png(const Framebuffer& fb) {
    const int w = fb.width(), h = fb.height();
    std::vector<uint8_t> raw;
    raw.reserve(static_cast<std::size_t>(h) * (1 + 3 * w));
    for (int y = 0; y < h; ++y) {
        raw.push_back(0);  // filter: none
        for (int x = 0; x < w; ++x) {
            const uint32_t p = fb.at(x, y);
            raw.push_back(static_cast<uint8_t>(p >> 16));
            raw.push_back(static_cast<uint8_t>(p >> 8));
            raw.push_back(static_cast<uint8_t>(p));
        }
    }
    uLongf zlen = compressBound(raw.size());
    std::vector<uint8_t> z(zlen);
    if (compress2(z.data(), &zlen, raw.data(), raw.size(), 6) != Z_OK) {
        throw std::runtime_error("png: zlib compression failed");
    }
    z.resize(zlen);

    std::vector<uint8_t> out(std::begin(kSignature), std::end(kSignature));
    std::vector<uint8_t> ihdr;
    put32(ihdr, static_cast<uint32_t>(w));
    put32(ihdr, static_cast<uint32_t>(h));
    ihdr.insert(ihdr.end(), {8, 2, 0, 0, 0});  // 8-bit, RGB, deflate, filter 0, no interlace
    put_chunk(out, "IHDR", ihdr);
    put_chunk(out, "IDAT", z);
    put_chunk(out, "IEND", {});
    return out;
}

void write_png(const Framebuffer& fb, const std::string& path) {
    const auto bytes = encode_png(fb);
    std::ofstream f(path, std::ios::binary);
    f.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (!f) throw std::runtime_error("png: cannot write " + path);
}

Framebuffer decode_png(const std::vector<uint8_t>& in) {
    if (in.size() < 8 || std::memcmp(in.data(), kSignature, 8) != 0) {
        throw std::runtime_error("png: bad signature");
    }
    int w = 0, h = 0;
    std::vector<uint8_t> z;
    for (std::size_t pos = 8; pos + 12 <= in.size();) {
        const uint32_t len = get32(&in[pos]);
        if (pos + 12 + len > in.size()) throw std::runtime_error("png: truncated chunk");
        const uint8_t* type = &in[pos + 4];
        const uint8_t* body = &in[pos + 8];
        if (std::memcmp(type, "IHDR", 4) == 0) {
            w = static_cast<int>(get32(body));
            h = static_cast<int>(get32(body + 4));
            if (body[8] != 8 || body[9] != 2 || body[12] != 0) {
                throw std::runtime_error("png: only 8-bit RGB, non-interlaced is supported");
            }
        } else if (std::memcmp(type, "IDAT", 4) == 0) {
            z.insert(z.end(), body, body + len);
        }
        pos += 12 + len;
    }
    const std::size_t stride = 1 + 3 * static_cast<std::size_t>(w);
    std::vector<uint8_t> raw(stride * h);
    uLongf rawlen = raw.size();
    if (w <= 0 || h <= 0 || uncompress(raw.data(), &rawlen, z.data(), z.size()) != Z_OK ||
        rawlen != raw.size()) {
        throw std::runtime_error("png: bad image data");
    }
    Framebuffer fb(w, h);
    for (int y = 0; y < h; ++y) {
        const uint8_t* r = &raw[y * stride];
        if (r[0] != 0) throw std::runtime_error("png: only filter type 0 is supported");
        for (int x = 0; x < w; ++x) {
            fb.row(y)[x] = (uint32_t{r[1 + 3 * x]} << 16) | (uint32_t{r[2 + 3 * x]} << 8) | r[3 + 3 * x];
        }
    }
    return fb;
}

Framebuffer read_png(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("png: cannot read " + path);
    return decode_png(std::vector<uint8_t>(std::istreambuf_iterator<char>(f), {}));
}

}  // namespace ssp
