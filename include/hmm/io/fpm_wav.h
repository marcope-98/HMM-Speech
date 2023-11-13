#ifndef FPM_WAV_H_
#define FPM_WAV_H_
#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <immintrin.h>

// TODO: proper error handling (improper use of assert)
// TODO: set preprocessor macro to enable sse
// TODO: cmake check compiler support for sse

// clang-format off
#define FPM_WAV_RIFF 0x46464952 // "RIFF" reversed for big endian
#define FPM_WAV_WAVE 0x45564157 // "WAVE" reversed for big endian
#define FPM_WAV_FMT  0x20746d66 // "fmt " reversed for big endian
#define FPM_WAV_DATA 0x61746164 // "data" reversed for big endian
#define FPM_WAV_FACT 0x74636166 // "fact" reversed for big endian

#define FPM_WAVE_FORMAT_PCM        0x0001
#define FPM_WAVE_FORMAT_IEEE_FLOAT 0x0003
#define FPM_WAVE_FORMAT_ALAW       0x0006
#define FPM_WAVE_FORMAT_MULAW      0x0007
#define FPM_WAVE_FORMAT_EXTENSIBLE 0xFFFE
// clang-format on

typedef struct
{
  uint32_t id   = 0;
  uint32_t size = 0;
} fpm_wav_chunk;

typedef struct
{
  uint8_t *     at    = NULL;
  uint8_t *     end   = NULL;
  fpm_wav_chunk chunk = {0};
} fpm_wav_iterator;

typedef struct
{
  uint32_t riff = 0;
  uint32_t size = 0;
  uint32_t wave = 0;
} fpm_riff_header;

typedef struct
{
  uint16_t wFormatTag          = 0;
  uint16_t nChannels           = 0;
  uint32_t nSamplesPerSec      = 0;
  uint32_t nAvgBytesPerSec     = 0;
  uint16_t nBlockAlign         = 0;
  uint16_t wBitsPerSample      = 0;
  uint16_t cbSize              = 0;
  uint16_t wValidBitsPerSample = 0;
  uint32_t dwChannelMask       = 0;
  uint8_t  subFormat[16]       = {0};
} fpm_fmt;

typedef struct
{
  uint32_t dwSampleLength = 0;
} fpm_fact;

typedef struct
{
  fpm_fmt  fmt      = {0};
  fpm_fact fact     = {0};
  float *  data     = NULL;
  size_t   nSamples = {0}; // number of samples per channel
} fpm_audio;

float *fpm_wav_load(const char *filename, size_t *samples, size_t *channels, size_t *hz);

fpm_wav_iterator fpm_wav_get_chunk(fpm_wav_iterator iter);
fpm_wav_iterator fpm_wav_get_chunk_header(fpm_wav_iterator it);
fpm_wav_iterator fpm_wav_init_iterator(const void *buffer, size_t size);

void             fpm_wav_parse(const void *buffer, size_t size, fpm_audio *audio);
fpm_wav_iterator fpm_wav_parse_header(fpm_wav_iterator it);
void             fpm_wav_parse_fmt__chunk(fpm_wav_iterator it, fpm_audio *audio);
void             fpm_wav_parse_fact_chunk(fpm_wav_iterator it, fpm_audio *audio);
void             fpm_wav_parse_data_chunk(fpm_wav_iterator it, fpm_audio *audio);

void fpm_wav_pcm(fpm_wav_iterator it, fpm_audio *audio);
void fpm_wav_pcm_8bits(fpm_wav_iterator it, fpm_audio *audio);
void fpm_wav_pcm_16bits(fpm_wav_iterator it, fpm_audio *audio);
void fpm_wav_pcm_24bits(fpm_wav_iterator it, fpm_audio *audio);
void fpm_wav_pcm_32bits(fpm_wav_iterator it, fpm_audio *audio);

void   fpm_wav_info(const char *filename);
void   fpm_wav_free(float *data);
size_t fpm_wav_align(size_t value, size_t bytes);
#endif // FPM_WAV_H_

#ifdef FPM_WAV_IMPLEMENTATION

void *fpm_wav_read_file(const char *filename, size_t *size)
{
  void *    buffer = NULL;
  FILE *    file   = fopen(filename, "rb");
  long      fsize;
  size_t    unused;
  fpm_audio audio;

  if (file == NULL) goto error;

  // get size of file and store content in buffer
  if (fseek(file, 0, SEEK_END) < 0) goto error;

  fsize = ftell(file);
  if (fsize < 0) goto error;

  buffer = malloc(fsize);
  if (buffer == NULL) goto error;

  rewind(file);

  unused = fread(buffer, 1, fsize, file);
  assert(unused == (size_t)fsize);
  if (ferror(file)) goto error;

  // close file
  fclose(file);

  *size = (size_t)fsize;
  return buffer;

error:
  if (file) fclose(file);
  if (buffer) free(buffer);
  return NULL;
}

float *fpm_wav_load(const char *filename, size_t *samples, size_t *channels, size_t *hz)
{
  size_t fsize  = 0;
  void * buffer = fpm_wav_read_file(filename, &fsize);
  if (buffer == NULL) return NULL;
  // parse
  fpm_audio audio = {0};
  fpm_wav_parse(buffer, fsize, &audio);
  *samples  = audio.nSamples;
  *channels = (size_t)audio.fmt.nChannels;
  *hz       = (size_t)audio.fmt.nSamplesPerSec;

  // free buffer
  free(buffer);

  return audio.data;
}

void fpm_wav_parse(const void *buffer, size_t size, fpm_audio *audio)
{
  // parse header
  fpm_wav_iterator it = fpm_wav_init_iterator(buffer, size);
  it                  = fpm_wav_parse_header(it);

  for (it = fpm_wav_get_chunk(it);
       it.at < it.end;
       it = fpm_wav_get_chunk(it))
  {
    it = fpm_wav_get_chunk_header(it);
    switch (it.chunk.id)
    {
      case FPM_WAV_FMT:
        fpm_wav_parse_fmt__chunk(it, audio);
        break;
      case FPM_WAV_DATA:
        fpm_wav_parse_data_chunk(it, audio);
        it.chunk.size += (it.chunk.size % 2); // padding byte if nSamples is odd
        break;
      case FPM_WAV_FACT:
        fpm_wav_parse_fact_chunk(it, audio);
        break;
      default:
        break;
    }
  }
  assert(it.at == it.end);
}

fpm_wav_iterator fpm_wav_init_iterator(const void *buffer, size_t size)
{
  fpm_wav_iterator it = {0};
  it.at               = (uint8_t *)(buffer);
  it.end              = (uint8_t *)(buffer) + size;
  return it;
}

fpm_wav_iterator fpm_wav_get_chunk(fpm_wav_iterator it)
{
  it.at += it.chunk.size;
  return it;
}

fpm_wav_iterator fpm_wav_get_chunk_header(fpm_wav_iterator it)
{
  it.chunk = *(fpm_wav_chunk *)it.at;
  it.at += sizeof(fpm_wav_chunk);
  return it;
}

fpm_wav_iterator fpm_wav_parse_header(fpm_wav_iterator it)
{
  fpm_riff_header header = *(fpm_riff_header *)(it.at);
  // TODO: these should not be asserts
  assert(header.riff == FPM_WAV_RIFF);
  assert(header.size == (uint32_t)(it.end - it.at) - 8u);
  assert(header.wave == FPM_WAV_WAVE);
  it.at += sizeof(fpm_riff_header);
  return it;
}

void fpm_wav_parse_fmt__chunk(fpm_wav_iterator it, fpm_audio *audio)
{
  _mm_storeu_si128((__m128i *)(&(audio->fmt)), _mm_loadu_si128((__m128i *)(it.at)));
  if (it.chunk.size == 40)
  {
    audio->fmt.cbSize              = *(uint16_t *)(it.at + 16);
    audio->fmt.wValidBitsPerSample = *(uint16_t *)(it.at + 17);
    audio->fmt.dwChannelMask       = *(uint32_t *)(it.at + 18);
    _mm_storeu_si128((__m128i *)audio->fmt.subFormat, _mm_loadu_si128((__m128i *)(it.at + 19)));
  }
}

void fpm_wav_parse_fact_chunk(fpm_wav_iterator it, fpm_audio *audio)
{
  audio->fact.dwSampleLength = *(uint32_t *)it.at;
}

void fpm_wav_parse_data_chunk(fpm_wav_iterator it, fpm_audio *audio)
{
  audio->nSamples = it.chunk.size / (audio->fmt.wBitsPerSample >> 3);
  audio->data     = (float *)malloc(fpm_wav_align(audio->nSamples * sizeof(float), 16));

  switch (audio->fmt.wFormatTag)
  {
    case FPM_WAVE_FORMAT_PCM:
      fpm_wav_pcm(it, audio);
      break;
    case FPM_WAVE_FORMAT_IEEE_FLOAT:
    case FPM_WAVE_FORMAT_ALAW:
    case FPM_WAVE_FORMAT_MULAW:
    case FPM_WAVE_FORMAT_EXTENSIBLE:
    default:
      printf("Unsupported audio format\n");
      fpm_wav_free(audio->data);
      audio->nSamples = 0;
  }
}

void fpm_wav_info(const char *filename)
{
  size_t fsize  = 0;
  void * buffer = fpm_wav_read_file(filename, &fsize);

  // parse
  fpm_audio audio = {0};
  fpm_wav_parse(buffer, fsize, &audio);

  switch (audio.fmt.wFormatTag)
  {
    case FPM_WAVE_FORMAT_PCM:
      printf("Audio format:       PCM\n");
      break;
    case FPM_WAVE_FORMAT_IEEE_FLOAT:
      printf("Audio format:       IEEE float\n");
      break;
    case FPM_WAVE_FORMAT_ALAW:
      printf("Audio format:       A-law\n");
      break;
    case FPM_WAVE_FORMAT_MULAW:
      printf("Audio format:       mu-law\n");
      break;
    case FPM_WAVE_FORMAT_EXTENSIBLE:
      printf("Audio format:       Extensible\n");
      break;
    default:
      printf("Audio format:       Unknown\n");
  }

  printf("Number of channels: %d\n", audio.fmt.nChannels);
  printf("Sample rate:        %d\n", audio.fmt.nSamplesPerSec);
  printf("Byte rate:          %d\n", audio.fmt.nAvgBytesPerSec);
  printf("Block align:        %d\n", audio.fmt.nBlockAlign);
  printf("Bits per sample:    %d\n", audio.fmt.wBitsPerSample);
}

void fpm_wav_free(float *data)
{
  free(data);
  data = NULL;
}

size_t fpm_wav_align(size_t value, size_t bytes)
{
  return (value + bytes - 1) & ~(bytes - 1);
}

void fpm_wav_pcm(fpm_wav_iterator it, fpm_audio *audio)
{
  switch (audio->fmt.wBitsPerSample)
  {
    case 8:
      fpm_wav_pcm_8bits(it, audio);
      break;
    case 16:
      fpm_wav_pcm_16bits(it, audio);
      break;
    case 24:
      fpm_wav_pcm_24bits(it, audio);
      break;
    case 32:
      fpm_wav_pcm_32bits(it, audio);
      break;
    default:
      printf("Unsupported Bits per sample\n");
  }
}

void fpm_wav_pcm_8bits(fpm_wav_iterator it, fpm_audio *audio)
{
#if 0
  __m128i  value_si128;
  __m128   value_ps;
  uint8_t *src    = (uint8_t *)it.at;
  float *  dst    = audio->data;
  __m128   sub    = _mm_set1_ps(128.f);
  __m128   factor = _mm_set1_ps(1.f / 127.5f);

  size_t optim = (audio->nSamples > 16 ? 0 : audio->nSamples - 16);

  for (size_t i = 0; i < optim; i += 16)
  {
    value_si128 = _mm_loadu_si128((__m128i *)(src + i));
    value_ps    = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(value_si128));
    _mm_store_ps(dst + i + 0, _mm_mul_ps(_mm_sub_ps(value_ps, sub), factor));
    value_ps = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_bsrli_si128(value_si128, 4)));
    _mm_store_ps(dst + i + 4, _mm_mul_ps(_mm_sub_ps(value_ps, sub), factor));
    value_ps = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_bsrli_si128(value_si128, 8)));
    _mm_store_ps(dst + i + 8, _mm_mul_ps(_mm_sub_ps(value_ps, sub), factor));
    value_ps = _mm_cvtepi32_ps(_mm_cvtepu8_epi32(_mm_bsrli_si128(value_si128, 12)));
    _mm_store_ps(dst + i + 12, _mm_mul_ps(_mm_sub_ps(value_ps, sub), factor));
  }

  for (size_t i = optim; i < audio->nSamples; ++i)
    audio->data[i] = (float(src[i]) - 128.f) / 127.5f;

#else
  uint8_t *src = (uint8_t *)it.at;
  for (size_t i = 0; i < audio->nSamples; ++i)
    audio->data[i] = (float(src[i]) - 128.f) / 127.5f;
#endif
}

void fpm_wav_pcm_16bits(fpm_wav_iterator it, fpm_audio *audio)
{
#if 0
  __m128i   value_si128;
  __m128    value_ps;
  uint16_t *src    = (uint16_t *)it.at;
  float *   dst    = audio->data;
  __m128    factor = _mm_set1_ps(1.f / 32768.f);

  size_t optim = (audio->nSamples > 8 ? 0 : audio->nSamples - 8);

  for (size_t i = 0; i < optim; i += 8)
  {
    value_si128 = _mm_loadu_si128((__m128i *)(src + i));
    value_ps    = _mm_cvtepi32_ps(_mm_cvtepi16_epi32(value_si128));
    _mm_store_ps(dst + i + 0, _mm_mul_ps(value_ps, factor));
    value_ps = _mm_cvtepi32_ps(_mm_cvtepi16_epi32(_mm_bsrli_si128(value_si128, 8)));
    _mm_store_ps(dst + i + 4, _mm_mul_ps(value_ps, factor));
  }

  for (size_t i = optim; i < audio->nSamples; ++i)
    audio->data[i] = float(src[i]) / 32768.f;

#else
  int16_t *src = (int16_t *)it.at;
  for (size_t i = 0; i < audio->nSamples; ++i)
    audio->data[i] = float(src[i]) / 32768.f;
#endif
}

void fpm_wav_pcm_24bits(fpm_wav_iterator it, fpm_audio *audio)
{
#if 0
  __m128   value_ps;
  uint8_t *src    = (uint8_t *)it.at;
  float *  dst    = audio->data;
  __m128   factor = _mm_set1_ps(1.f / 2147483648.f);
  __m128i  mask   = _mm_set_epi32(0xff110109, 0xff080706, 0xff050403, 0xff020100);

  size_t optim = (audio->nSamples > 4 ? 0 : audio->nSamples - 4);

  for (size_t i = 0; i < optim; i += 4)
  {
    value_ps = _mm_cvtepi32_ps(_mm_shuffle_epi8(_mm_loadu_si128((__m128i *)src), mask));
    _mm_store_ps(dst + i, _mm_mul_ps(value_ps, factor));
    src += 12;
  }

  for (size_t i = optim; i < audio->nSamples; ++i)
  {
    int32_t temp   = *(int32_t *)(src + (i * 3)) & 0xffffff00;
    audio->data[i] = float(temp) / 2147483648.f;
  }

#else
  uint8_t *src = (uint8_t *)it.at;
  for (size_t i = 0; i < audio->nSamples; ++i)
  {
    int32_t temp   = *(int32_t *)(src + (i * 3)) & 0xffffff00;
    audio->data[i] = float(temp) / 2147483648.f;
  }
#endif
}

void fpm_wav_pcm_32bits(fpm_wav_iterator it, fpm_audio *audio)
{
#if 0
  int32_t *src    = (int32_t *)it.at;
  float *  dst    = audio->data;
  __m128   factor = _mm_set1_ps(1.f / 2147483648.f);

  size_t optim = (audio->nSamples > 4 ? 0 : audio->nSamples - 4);

  for (size_t i = 0; i < optim; i += 4)
    _mm_store_ps(dst + i, _mm_mul_ps(_mm_cvtepi32_ps(_mm_loadu_si128((__m128i *)(src + i))), factor));
  for (size_t i = optim; i < audio->nSamples; ++i)
    audio->data[i] = float(src[i]) / 2147483648.f;
#else
  int32_t *src = (int32_t *)it.at;
  for (size_t i = 0; i < audio->nSamples; ++i)
    audio->data[i] = float(src[i]) / 2147483648.f;
#endif
}

#endif // FPM_WAV_IMPLEMENTATION