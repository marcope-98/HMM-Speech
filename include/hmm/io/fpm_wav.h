#ifndef FPM_WAV_H_
#define FPM_WAV_H_
#include <assert.h>
#include <errno.h>
#include <immintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FPM_WAV_RIFF 0x46464952 // reversed for big endian
#define FPM_WAV_WAVE 0x45564157 // reversed for big endian
#define FPM_WAV_FMT 0x20746d66  // reversed for big endian
#define FPM_WAV_DATA 0x61746164 // reversed for big endian
#define FPM_WAV_FACT 0x74636166 // reversed for big endian

// NOTE: These are here just to help extend library later, as of now only PCM is supported
#define FPM_WAVE_FORMAT_PCM 0x0001
#define FPM_WAVE_FORMAT_IEEE_FLOAT 0x0003
#define FPM_WAVE_FORMAT_ALAW 0x0006
#define FPM_WAVE_FORMAT_MULAW 0x0007
#define FPM_WAVE_FORMAT_EXTENSIBLE 0xFFFE

// http://soundfile.sapp.org/doc/WaveFormat/

struct fpm_fmt
{
  unsigned short wFormatTag          = {0};
  unsigned short nChannels           = {0};
  unsigned int   nSamplesPerSec      = {0};
  unsigned int   nAvgBytesPerSec     = {0};
  unsigned short nBlockAlign         = {0};
  unsigned short wBitsPerSample      = {0};
  unsigned short cbSize              = {0};
  unsigned short wValidBitsPerSample = {0};
  unsigned int   dwChannelMask       = {0};
  unsigned char  SubFormat[16]       = {0};
};

struct fpm_fact
{
  unsigned int dwSampleLength = {0};
};

struct fpm_audio
{
  fpm_fmt  fmt    = {0};
  fpm_fact fact   = {0};
  float *  data   = NULL;
  size_t   stride = {0}; // number of samples per channel
};

size_t fpm_get_file_size(FILE *file);

int fpm_wav_load(const char *filename, fpm_audio *audio);

void   fpm_wav_parse(const unsigned int *buffer, size_t size, fpm_audio *audio);
size_t fpm_wav_parse_header(const unsigned int *buffer, size_t size);
size_t fpm_wav_parse_fmt__chunk(const unsigned int *buffer, fpm_audio *audio);
size_t fpm_wav_parse_data_chunk(const unsigned int *buffer, fpm_audio *audio);
size_t fpm_wav_parse_fact_chunk(const unsigned int *buffer, fpm_audio *audio);
size_t fpm_wav_skip_chunk(const unsigned int *buffer);

void fpm_wave_pcm(const unsigned int *buffer, size_t size, fpm_audio *audio);
void fpm_wave_pcm_8bits(const unsigned int *buffer, size_t size, fpm_audio *audio);
void fpm_wave_pcm_16bits(const unsigned int *buffer, size_t size, fpm_audio *audio);
void fpm_wave_pcm_24bits(const unsigned int *buffer, size_t size, fpm_audio *audio);
void fpm_wave_pcm_32bits(const unsigned int *buffer, size_t size, fpm_audio *audio);
void fpm_wave_ieee_float(const unsigned int *buffer, size_t size, fpm_audio *audio);
void fpm_wave_alaw(const unsigned int *buffer, size_t size, fpm_audio *audio);
void fpm_wave_mulaw(const unsigned int *buffer, size_t size, fpm_audio *audio);
void fpm_wave_extensible(const unsigned int *buffer, size_t size, fpm_audio *audio);

void fpm_wav_info(fpm_audio *audio);
void fpm_wav_free(fpm_audio *audio);
#endif // FPM_WAV_H_

#ifdef FPM_WAV_IMPLEMENTATION

size_t fpm_get_file_size(FILE *file)
{
  fseek(file, 0L, SEEK_END);
  size_t size = ftell(file);
  rewind(file);
  return size;
}

int fpm_wav_load(const char *filename, fpm_audio *audio)
{
  // open file
  unsigned int *buffer = NULL;
  FILE *        file   = fopen(filename, "rb");
  long          fsize;
  size_t        unused;

  if (file == NULL)
    goto error;

  // get size of file and store content in buffer buffer
  if (fseek(file, 0, SEEK_END) < 0)
    goto error;

  fsize = ftell(file);
  if (fsize < 0)
    goto error;

  buffer = (unsigned int *)malloc(fsize);
  if (buffer == NULL)
    goto error;

  rewind(file);

  unused = fread(buffer, 1, fsize, file);
  assert(unused == (size_t)fsize);
  if (ferror(file))
    goto error;

  // close file
  fclose(file);

  // parse
  fpm_wav_parse(buffer, fsize, audio);

  // free buffer
  free(buffer);

  return 0;

error:
  if (file) fclose(file);
  if (buffer) free(buffer);
  return 1;
}

void fpm_wav_parse(const unsigned int *buffer, size_t size, fpm_audio *audio)
{
  size_t offset = 0;
  offset += fpm_wav_parse_header(buffer + offset, size);

  while (offset < (size >> 2))
  {
    switch (*(buffer + (offset++)))
    {
      case FPM_WAV_FMT: // fmt  chunk
        offset += fpm_wav_parse_fmt__chunk(buffer + offset, audio);
        break;
      case FPM_WAV_DATA: // data chunk
        offset += fpm_wav_parse_data_chunk(buffer + offset, audio);
        break;
      case FPM_WAV_FACT: // fact chunk
        offset += fpm_wav_parse_fact_chunk(buffer + offset, audio);
        break;
      default: // skip
        offset += fpm_wav_skip_chunk(buffer + offset);
    }
  }

  assert(offset == (size >> 2) && "Size does not match: data subchunk might contain additional subchunks");
}

size_t fpm_wav_parse_header(const unsigned int *buffer, size_t size)
{
  const __m128i header = _mm_set_epi32(FPM_WAV_FMT, FPM_WAV_WAVE, size - 8, FPM_WAV_RIFF);

  // _mm_test_all_zeros returns 1 if they are equal
  __m128i testeq = _mm_xor_si128(header, _mm_loadu_si128((__m128i *)buffer));
  if (_mm_test_all_zeros(testeq, testeq) == 0)
  {
    fprintf(stderr, "[ERROR] Wav file header is corrupted\n");
    exit(1);
  }
  return 3;
}

size_t fpm_wav_parse_fmt__chunk(const unsigned int *buffer, fpm_audio *audio)
{
  // chunk size
  size_t subchunksize = *(buffer++);
  _mm_storeu_si128((__m128i *)(&(audio->fmt)), _mm_loadu_si128((__m128i *)buffer));
  buffer += 4;

  if (subchunksize == 40)
  {
    size_t mask                    = *(buffer++);
    audio->fmt.cbSize              = mask >> 4;
    audio->fmt.wValidBitsPerSample = mask & 0xFFFF;
    audio->fmt.dwChannelMask       = *(buffer++);
    _mm_storeu_si128((__m128i *)audio->fmt.SubFormat, _mm_loadu_si128((__m128i *)buffer));
  }

  return (1 + (subchunksize >> 2)); // skip extra parameters
}

size_t fpm_wav_parse_data_chunk(const unsigned int *buffer, fpm_audio *audio)
{
  // subchunk2size
  size_t subchunksize = *(buffer++);

  // the data depends on wFormatTag, channel number and wBitspersample
  audio->stride = subchunksize / (audio->fmt.nChannels * (audio->fmt.wBitsPerSample >> 3)); // number of samples
  audio->data   = (float *)malloc(audio->stride * audio->fmt.nChannels * sizeof(float));

  // clang-format off
  switch (audio->fmt.wFormatTag)
  {
    case FPM_WAVE_FORMAT_PCM:
      fpm_wave_pcm(buffer, subchunksize, audio); break;
    case FPM_WAVE_FORMAT_IEEE_FLOAT:
      fpm_wave_ieee_float(buffer, subchunksize, audio); break;
    case FPM_WAVE_FORMAT_ALAW:
      fpm_wave_alaw(buffer, subchunksize, audio); break;
    case FPM_WAVE_FORMAT_MULAW:
      fpm_wave_mulaw(buffer, subchunksize, audio); break;
    case FPM_WAVE_FORMAT_EXTENSIBLE:
      fpm_wave_extensible(buffer, subchunksize, audio); break;
    default:
      printf("Unsupported audio format\n");
  }
  // clang-format on

  return 1 + ((subchunksize + (subchunksize % 2)) >> 2); // consider pad byte
}

size_t fpm_wav_parse_fact_chunk(const unsigned int *buffer, fpm_audio *audio)
{
  // chunk size
  size_t subchunksize = *(buffer++);
  // dwSampleLength
  audio->fact.dwSampleLength = *buffer;
  return 1 + (subchunksize >> 2);
}

size_t fpm_wav_skip_chunk(const unsigned int *buffer) { return 1 + ((*buffer) >> 2); }

void fpm_wave_pcm(const unsigned int *buffer, size_t size, fpm_audio *audio)
{
  switch (audio->fmt.wBitsPerSample)
  {
    case 8:
      fpm_wave_pcm_8bits(buffer, size, audio);
      break;
    case 16:
      fpm_wave_pcm_16bits(buffer, size, audio);
      break;
    case 24:
      fpm_wave_pcm_24bits(buffer, size, audio);
      break;
    case 32:
      fpm_wave_pcm_32bits(buffer, size, audio);
      break;
    default:
      printf("Unsupported Bits per sample\n");
      return;
  }
  // if (audio->fmt.nChannels == 1 && audio->fmt.wBitsPerSample == 16)
  // {
  // short *src      = (short *)buffer;
  // size_t elements = size / (audio->fmt.wBitsPerSample >> 3);
  //
  // for (size_t i = 0; i < elements; ++i)
  // audio->data[i] = float(src[i]);
  // }
}
void fpm_wave_pcm_8bits(const unsigned int *buffer, size_t size, fpm_audio *audio)
{
  unsigned char *src      = (unsigned char *)buffer;
  size_t         elements = size;
  for (size_t i = 0; i < elements; ++i)
    audio->data[i] = ((float(src[i]) - 128.f) / 255.f) * 2.f; // centered in 0x80 and with 1 of amplitude
}
void fpm_wave_pcm_16bits(const unsigned int *buffer, size_t size, fpm_audio *audio)
{
  short *src      = (short *)buffer;
  size_t elements = size / (audio->fmt.wBitsPerSample >> 3);
  for (size_t i = 0; i < elements; ++i)
    audio->data[i] = float(src[i]) / float(32768.f);
}

void fpm_wave_pcm_24bits(const unsigned int *buffer, size_t size, fpm_audio *audio)
{
  int *  src      = (int *)buffer;
  size_t elements = size / (audio->fmt.wBitsPerSample >> 3);
  int    tmp;

  // elements must be divisible by 3
  assert(elements % 3 == 0 && "report as bug: function name fpm_wave_pcm_24bits");
  for (size_t i = 0, j = 0; i < elements; i += 3, j += 4)
  {
    // we process 3 elements at a time:
    // 0xaa aa aa bb
    // 0xbb bb cc cc
    // 0xcc dd dd dd

    tmp            = src[i];                                        // 0xaa aa aa bb
    audio->data[j] = float(tmp & 0xFFFFFF00) / float(2147483648.f); // 0xaa aa aa 00
    tmp <<= 24;                                                     // 0xbb 00 00 00

    tmp |= (src[i + 1] >> 8);                                           // 0xbb bb bb cc
    audio->data[j + 1] = float(tmp & 0xFFFFFF00) / float(2147483648.f); // 0xbb bb bb 00
    tmp                = src[i + 1] << 16;                              // 0xcc cc 00 00

    tmp |= (src[i + 2] >> 16);                                          // 0xcc cc cc dd
    audio->data[j + 2] = float(tmp & 0xFFFFFF00) / float(2147483648.f); // 0xcc cc cc 00
    tmp                = src[i + 2] << 8;                               // 0xdd dd dd 00

    audio->data[j + 3] = float(tmp & 0xFFFFFF00) / float(2147483648.f);
  }
}

void fpm_wave_pcm_32bits(const unsigned int *buffer, size_t size, fpm_audio *audio)
{
  int *  src      = (int *)buffer;
  size_t elements = size / (audio->fmt.wBitsPerSample >> 3);
  for (size_t i = 0; i < elements; ++i)
    audio->data[i] = float(src[i]) / float(2147483648.f);
}

void fpm_wave_ieee_float(const unsigned int *buffer, size_t size, fpm_audio *audio)
{
  (void)buffer;
  (void)size;
  (void)audio;
}
void fpm_wave_alaw(const unsigned int *buffer, size_t size, fpm_audio *audio)
{
  (void)buffer;
  (void)size;
  (void)audio;
}
void fpm_wave_mulaw(const unsigned int *buffer, size_t size, fpm_audio *audio)
{
  (void)buffer;
  (void)size;
  (void)audio;
}
void fpm_wave_extensible(const unsigned int *buffer, size_t size, fpm_audio *audio)
{
  (void)buffer;
  (void)size;
  (void)audio;
}

void fpm_wav_free(fpm_audio *audio) { free(audio->data); }

void fpm_wav_info(fpm_audio *audio)
{
  switch (audio->fmt.wFormatTag)
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

  printf("Number of channels: %d\n", audio->fmt.nChannels);
  printf("Sample rate:        %d\n", audio->fmt.nSamplesPerSec);
  printf("Byte rate:          %d\n", audio->fmt.nAvgBytesPerSec);
  printf("Block align:        %d\n", audio->fmt.nBlockAlign);
  printf("Bits per sample:    %d\n", audio->fmt.wBitsPerSample);
}

#endif