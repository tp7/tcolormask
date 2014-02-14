#include <vector>
#include <regex>
#include <future>
#include <stdint.h>
#include <Windows.h>
#pragma warning(disable: 4512 4244 4100)
#include <avisynth.h>
#pragma warning(default: 4512 4244 4100)
#include <xmmintrin.h>

using namespace std;

struct YUVPixel {
    uint8_t Y;
    uint8_t U;
    uint8_t V;

    uint32_t vector_y;
    uint32_t vector_u;
    uint32_t vector_v;
};


int round(float d) {
    return static_cast<int>(d + 0.5f);
}

template<int subsamplingX, int subsamplingY>
void processLut(uint8_t *pDstY, const uint8_t *pSrcY, const uint8_t *pSrcV, const uint8_t *pSrcU, int dstPitchY, int srcPitchY, int srcPitchUV, int width, int height, uint8_t *lutY, uint8_t *lutU, uint8_t *lutV) {
    for(int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            pDstY[x] = lutY[pSrcY[x]] & lutU[pSrcU[x/subsamplingX]] & lutV[pSrcV[x/subsamplingX]];
        }
        pSrcY += srcPitchY;
        if (y % subsamplingY == (subsamplingY-1)) {
            pSrcU += srcPitchUV;
            pSrcV += srcPitchUV;
        }
        pDstY += dstPitchY;
    }
}

template<int subsamplingX, int subsamplingY>
void processSse2(uint8_t *pDstY, const uint8_t *pSrcY, const uint8_t *pSrcV, const uint8_t *pSrcU, int dstPitchY, 
                 int srcPitchY, int srcPitchUV, int width, int height, const vector<YUVPixel>& colors, int vectorTolerance, int halfVectorTolerance) {
    for(int y = 0; y < height; ++y) {
        for (int x = 0; x < width; x+=16) {
            auto result_y = _mm_setzero_si128();
            auto result_u = _mm_setzero_si128();
            auto result_v = _mm_setzero_si128();

            auto srcY_v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(pSrcY+x));
            __m128i srcU_v, srcV_v;
            if (subsamplingX == 2) {
                srcU_v = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(pSrcU+x/subsamplingX));
                srcU_v = _mm_unpacklo_epi8(srcU_v, srcU_v);
                srcV_v = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(pSrcV+x/subsamplingX));
                srcV_v = _mm_unpacklo_epi8(srcV_v, srcV_v);
            } else {
                srcU_v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(pSrcU+x));
                srcV_v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(pSrcV+x));
            }

            for(auto &color: colors) {
                auto colorVector_y = _mm_set1_epi32(color.vector_y);
                auto colorVector_u = _mm_set1_epi32(color.vector_u);
                auto colorVector_v = _mm_set1_epi32(color.vector_v);
                /* absolute difference */
                auto maximum_y = _mm_max_epu8(srcY_v, colorVector_y);
                auto maximum_u = _mm_max_epu8(srcU_v, colorVector_u);
                auto maximum_v = _mm_max_epu8(srcV_v, colorVector_v);

                auto minimum_y = _mm_min_epu8(srcY_v, colorVector_y);
                auto minimum_u = _mm_min_epu8(srcU_v, colorVector_u);
                auto minimum_v = _mm_min_epu8(srcV_v, colorVector_v);

                auto diff_y = _mm_subs_epu8(maximum_y, minimum_y);
                auto diff_u = _mm_subs_epu8(maximum_u, minimum_u);
                auto diff_v = _mm_subs_epu8(maximum_v, minimum_v);
                /* comparing to tolerance */
                auto diff_tolerance_min_y = _mm_max_epu8(diff_y, _mm_set1_epi32(vectorTolerance));
                auto diff_tolerance_min_u = _mm_max_epu8(diff_u, _mm_set1_epi32(halfVectorTolerance));
                auto diff_tolerance_min_v = _mm_max_epu8(diff_v, _mm_set1_epi32(halfVectorTolerance));

                auto passed_y = _mm_cmpeq_epi8(diff_y, diff_tolerance_min_y);
                auto passed_u = _mm_cmpeq_epi8(diff_u, diff_tolerance_min_u);
                auto passed_v = _mm_cmpeq_epi8(diff_v, diff_tolerance_min_v);
                /* inverting to get "lower" instead of "lower or equal" */
                passed_y = _mm_andnot_si128(passed_y, _mm_set1_epi32(0xFFFFFFFF));
                passed_u = _mm_andnot_si128(passed_u, _mm_set1_epi32(0xFFFFFFFF));
                passed_v = _mm_andnot_si128(passed_v, _mm_set1_epi32(0xFFFFFFFF));

                result_y = _mm_or_si128(result_y, passed_y);
                result_u = _mm_or_si128(result_u, passed_u);
                result_v = _mm_or_si128(result_v, passed_v);
            }
            result_y = _mm_and_si128(result_y, result_u);
            result_y = _mm_and_si128(result_y, result_v);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(pDstY+x), result_y);
        }
        pSrcY += srcPitchY;
        if (y % subsamplingY == (subsamplingY-1)) {
            pSrcU += srcPitchUV;
            pSrcV += srcPitchUV;
        }
        pDstY += dstPitchY;
    }
}


auto lutYv12 = &processLut<2, 2>;
auto lutYv24 = &processLut<1, 1>;
auto lutYv16 = &processLut<2, 1>;
auto sse2Yv12 = &processSse2<2, 2>;
auto sse2Yv24 = &processSse2<1, 1>;
auto sse2Yv16 = &processSse2<2, 1>;

class TColorMask : public GenericVideoFilter {
public:
    TColorMask(PClip child, vector<int> colors, int tolerance, bool bt601, bool grayscale, int lutthr, bool mt, IScriptEnvironment* env);
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
    
    int __stdcall SetCacheHints(int cachehints, int frame_range) override {
        return cachehints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0;
    }

private:
    void buildLuts();

    void process(uint8_t *dstY_ptr, const uint8_t *srcY_ptr, const uint8_t *srcV_ptr, const uint8_t *srcU_ptr, int dst_pitch_y, int src_pitch_y, int src_pitch_uv, int width, int height);

    vector<YUVPixel> colors_;
    int tolerance_;
    bool grayscale_;
    bool mt_;
    int prefer_lut_thresh_;
    int subsamplingY_;
    uint32_t vector_tolerance_;
    uint32_t vector_half_tolerance_;

    uint8_t lut_y[256];
    uint8_t lut_u[256];
    uint8_t lut_v[256];
    
    decltype(lutYv12) lutFunction_;
    decltype(sse2Yv12) sse2Function_;
};

TColorMask::TColorMask(PClip child, vector<int> colors, int tolerance, bool bt601, bool grayscale, int lutthr, bool mt, IScriptEnvironment* env) 
    : GenericVideoFilter(child), tolerance_(tolerance), grayscale_(grayscale), prefer_lut_thresh_(lutthr), mt_(mt), vector_tolerance_(0), vector_half_tolerance_(0) {
    if (vi.IsYV24()) {
        subsamplingY_ = 1;
        lutFunction_ = lutYv24;
        sse2Function_ = sse2Yv24;
    } else if (vi.IsYV12()) {
        subsamplingY_ = 2;
        lutFunction_ = lutYv12;
        sse2Function_ = sse2Yv12;
    } else if (vi.IsYV16()) {
        subsamplingY_ = 1;
        lutFunction_ = lutYv16;
        sse2Function_ = sse2Yv16;
    } else {
        env->ThrowError("Only YV12, YV16 and YV24 are supported!");
    }
    float kR = bt601 ? 0.299f : 0.2126f;
    float kB = bt601 ? 0.114f : 0.0722f;

    for (auto color: colors) {
        YUVPixel p;
        memset(&p, 0, sizeof(p));
        float r = static_cast<float>((color & 0xFF0000) >> 16) / 255.0f;
        float g = static_cast<float>((color & 0xFF00) >> 8) / 255.0f;
        float b = static_cast<float>(color & 0xFF) / 255.0f;

        float y = kR*r + (1-kR-kB)*g + kB*b;
        p.U = 128 + round(112.0f*(b-y)/(1-kB));
        p.V = 128 + round(112.0f*(r-y)/(1-kR));
        p.Y = 16 + round(219.0f*y);
        for (int i = 0; i < 4; i++) {
            p.vector_y |= p.Y << (8*i);
            p.vector_u |= p.U << (8*i);
            p.vector_v |= p.V << (8*i);
        }

        colors_.push_back(p);
    }

    for (int i = 0; i < 4; i++) {
        vector_tolerance_ |=  tolerance << (8*i);
        vector_half_tolerance_ |= vector_half_tolerance_ | (tolerance / 2) << (8*i);
    }
    
    if (((child->GetVideoInfo().width % 16) != 0) || (colors_.size() > prefer_lut_thresh_)) {
        buildLuts();
    }
}

void TColorMask::buildLuts() {
    for (int i = 0; i < 256; ++i) {
        uint8_t val_y = 0;
        uint8_t val_u = 0;
        uint8_t val_v = 0;
        for (auto &color: colors_) {
            val_y |= ((abs(i - color.Y) < tolerance_) ? 255 : 0);
            val_u |= ((abs(i - color.U) < (tolerance_ / 2)) ? 255 : 0);
            val_v |= ((abs(i - color.V) < (tolerance_ / 2)) ? 255 : 0);
        }
        lut_y[i] = val_y;
        lut_u[i] = val_u;
        lut_v[i] = val_v;
    }
}

PVideoFrame TColorMask::GetFrame(int n, IScriptEnvironment* env) {
   PVideoFrame src = child->GetFrame(n,env);
   auto dst = env->NewVideoFrame(child->GetVideoInfo());

   int width = src->GetRowSize(PLANAR_Y);
   int height = src->GetHeight(PLANAR_Y);

   const uint8_t *srcY_ptr = src->GetReadPtr(PLANAR_Y);
   const uint8_t *srcU_ptr = src->GetReadPtr(PLANAR_U);
   const uint8_t *srcV_ptr = src->GetReadPtr(PLANAR_V);
   int src_pitch_y = src->GetPitch(PLANAR_Y);
   int src_pitch_uv = src->GetPitch(PLANAR_U);
   
   uint8_t *dstY_ptr = dst->GetWritePtr(PLANAR_Y);
   int dst_pitch_y = dst->GetPitch(PLANAR_Y);
   
   if (mt_) {
       //async seems to be threadpool'ed on windows, creating threads is less efficient
       auto thread2 = std::async(launch::async, [=] { 
           process(dstY_ptr, 
               srcY_ptr, 
               srcV_ptr, 
               srcU_ptr, 
               dst_pitch_y, 
               src_pitch_y,
               src_pitch_uv,
               width,
               height/2); 
       });
       process(dstY_ptr + (dst_pitch_y*height/2), 
           srcY_ptr + (src_pitch_y*height/2), 
           srcV_ptr + (src_pitch_uv*height/(2*subsamplingY_)), 
           srcU_ptr + (src_pitch_uv*height/(2*subsamplingY_)), 
           dst_pitch_y, 
           src_pitch_y, 
           src_pitch_uv, 
           width, 
           height/2); 
       thread2.wait();
   } else {
       process(dstY_ptr, srcY_ptr, srcV_ptr, srcU_ptr, dst_pitch_y, src_pitch_y, src_pitch_uv, width, height);
   }

   if (grayscale_) {
       memset(dst->GetWritePtr(PLANAR_U), 128, dst->GetPitch(PLANAR_U) * dst->GetHeight(PLANAR_U));
       memset(dst->GetWritePtr(PLANAR_V), 128, dst->GetPitch(PLANAR_V) * dst->GetHeight(PLANAR_V));
   }
   return dst;
}

void TColorMask::process(uint8_t *dstY_ptr, const uint8_t *srcY_ptr, const uint8_t *srcV_ptr, const uint8_t *srcU_ptr, int dst_pitch_y, int src_pitch_y, int src_pitch_uv, int width, int height) {
    if (colors_.size() > prefer_lut_thresh_) {
        lutFunction_(dstY_ptr, srcY_ptr, srcV_ptr, srcU_ptr, dst_pitch_y, src_pitch_y, src_pitch_uv, width, height, lut_y, lut_u, lut_v);
        return;
    }
    int border = width % 16;

    sse2Function_(dstY_ptr, srcY_ptr, srcV_ptr, srcU_ptr, dst_pitch_y, src_pitch_y, src_pitch_uv, width - border , height, colors_, vector_tolerance_, vector_half_tolerance_);
    if (border != 0) {
        lutFunction_(dstY_ptr + width - border, srcY_ptr + width - border, srcV_ptr + width - border, srcU_ptr + width - border, dst_pitch_y, src_pitch_y, src_pitch_uv, border, height, lut_y, lut_u, lut_v);
    }
}

int avisynthStringToInt(const string &str) {
    if (str[0] == '$') {
        auto substr = str.substr(1, str.length());
        return strtol(substr.c_str(), 0, 16);
    }
    return strtol(str.c_str(), 0, 10);
}

AVSValue __cdecl CreateTColorMask(AVSValue args, void*, IScriptEnvironment* env) 
{
    enum { CLIP, COLORS, TOLERANCE, BT601, GRAYSCALE, LUTTHR, MT };

    string str = args[COLORS].AsString("");
    
    std::regex e ("(/\\*.*\\*/|//.*$)");   //comments
    str = regex_replace(str, e, "");

    vector<int> colors;
    regex rgx("(\\$?[\\da-fA-F]+)");

    sregex_token_iterator iter(str.cbegin(), str.cend(), rgx, 1);
    sregex_token_iterator end;
    for( ; iter != end; ++iter ) {
        try {
            colors.push_back(avisynthStringToInt(iter->str()));
        } catch(...) {
            env->ThrowError("Parsing error");
        }
    }
    
    return new TColorMask(args[CLIP].AsClip(), colors, args[TOLERANCE].AsInt(10), args[BT601].AsBool(false), 
        args[GRAYSCALE].AsBool(false), args[LUTTHR].AsInt(9), args[MT].AsBool(false), env);
}

const AVS_Linkage *AVS_linkage = nullptr;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors) {
    AVS_linkage = vectors;

    env->AddFunction("tcolormask", "c[colors]s[tolerance]i[bt601]b[gray]b[lutthr]i[mt]b", CreateTColorMask, 0);
    return "Why are you looking at this?";
}
