#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#define INBUFSIZ   1024
#define OUTBUFSIZ  1024

z_stream z;

char inbuf[INBUFSIZ];
char outbuf[OUTBUFSIZ];

void write_zfile(FILE *fin, FILE *fout)
{
    int count, flush, status;

    /* すべてのメモリ管理をライブラリに任せる */
    z.zalloc = Z_NULL;
    z.zfree = Z_NULL;
    z.opaque = Z_NULL;

    /* 初期化 */
    /* 第2引数は圧縮の度合。0〜9 の範囲の整数で，0 は無圧縮 */
    /* Z_DEFAULT_COMPRESSION (= 6) が標準 */
    if (deflateInit(&z, Z_DEFAULT_COMPRESSION) != Z_OK) {
        fprintf(stderr, "deflateInit: %s\n", (z.msg) ? z.msg : "???");
        exit(1);
    }

    z.avail_in = 0;             /* 入力バッファ中のデータのバイト数 */
    z.next_out = outbuf;        /* 出力ポインタ */
    z.avail_out = OUTBUFSIZ;    /* 出力バッファのサイズ */

    /* 通常は deflate() の第2引数は Z_NO_FLUSH にして呼び出す */
    flush = Z_NO_FLUSH;

    while (1) {
        if (z.avail_in == 0) {  /* 入力が尽きれば */
            z.next_in = inbuf;  /* 入力ポインタを入力バッファの先頭に */
            z.avail_in = fread(inbuf, 1, INBUFSIZ, fin); /* データを読み込む */

            /* 入力が最後になったら deflate() の第2引数は Z_FINISH にする */
            if (z.avail_in < INBUFSIZ) {
				flush = Z_FINISH;
			}
        }
        status = deflate(&z, flush); /* 圧縮する */
		/*
		printf( "next_in = %x, avail_in = %d, total_in = %d, avail_out = %d, out = %d.\n", z.next_in, z.avail_in, z.total_in, z.avail_out, z.total_out );
					         
							 */
		
        if (status == Z_STREAM_END) break; /* 完了 */
        if (status != Z_OK) {   /* エラー */
            fprintf(stderr, "deflate: %s\n", (z.msg) ? z.msg : "???");
            exit(1);
        }
        if (z.avail_out == 0) { /* 出力バッファが尽きれば */
            /* まとめて書き出す */
            if (fwrite(outbuf, 1, OUTBUFSIZ, fout) != OUTBUFSIZ) {
                fprintf(stderr, "Write error\n");
                exit(1);
            }
            z.next_out = outbuf; /* 出力バッファ残量を元に戻す */
            z.avail_out = OUTBUFSIZ; /* 出力ポインタを元に戻す */
        }
    }

    /* 残りを吐き出す */
    if ((count = OUTBUFSIZ - z.avail_out) != 0) {
        if (fwrite(outbuf, 1, count, fout) != count) {
            fprintf(stderr, "Write error\n");
            exit(1);
        }
    }

    /* 後始末 */
    if (deflateEnd(&z) != Z_OK) {
        fprintf(stderr, "deflateEnd: %s\n", (z.msg) ? z.msg : "???");
        exit(1);
    }
}

void read_zfile(FILE *fin, FILE *fout)
{
    int count, status;

    /* すべてのメモリ管理をライブラリに任せる */
    z.zalloc = Z_NULL;
    z.zfree = Z_NULL;
    z.opaque = Z_NULL;

    /* 初期化 */
    z.next_in = Z_NULL;
    z.avail_in = 0;
    if (inflateInit(&z) != Z_OK) {
        fprintf(stderr, "inflateInit: %s\n", (z.msg) ? z.msg : "???");
        exit(1);
    }

    z.next_out = outbuf;        /* 出力ポインタ */
    z.avail_out = OUTBUFSIZ;    /* 出力バッファ残量 */
    status = Z_OK;

    while (status != Z_STREAM_END) {
        if (z.avail_in == 0) {  /* 入力残量がゼロになれば */
            z.next_in = inbuf;  /* 入力ポインタを元に戻す */
            z.avail_in = fread(inbuf, 1, INBUFSIZ, fin); /* データを読む */
        }
        status = inflate(&z, Z_NO_FLUSH); /* 展開 */
		/*
		printf( "next_in = %x, avail_in = %d, total_in = %d, avail_out = %d, out = %d.\n", z.next_in, z.avail_in, z.total_in, z.avail_out, z.total_out );
		*/
        if (status == Z_STREAM_END) break; /* 完了 */
        if (status != Z_OK) {   /* エラー */
            fprintf(stderr, "inflate: %s\n", (z.msg) ? z.msg : "???");
            exit(1);
        }
        if (z.avail_out == 0) { /* 出力バッファが尽きれば */
            /* まとめて書き出す */
            if (fwrite(outbuf, 1, OUTBUFSIZ, fout) != OUTBUFSIZ) {
                fprintf(stderr, "Write error\n");
                exit(1);
            }
            z.next_out = outbuf; /* 出力ポインタを元に戻す */
            z.avail_out = OUTBUFSIZ; /* 出力バッファ残量を元に戻す */
        }
    }

    /* 残りを吐き出す */
    if ((count = OUTBUFSIZ - z.avail_out) != 0) {
        if (fwrite(outbuf, 1, count, fout) != count) {
            fprintf(stderr, "Write error\n");
            exit(1);
        }
    }

    /* 後始末 */
    if (inflateEnd(&z) != Z_OK) {
        fprintf(stderr, "inflateEnd: %s\n", (z.msg) ? z.msg : "???");
        exit(1);
    }
}

