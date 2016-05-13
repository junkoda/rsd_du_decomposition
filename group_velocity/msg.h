#ifndef MSG_H
#define MSG_H 1

#ifdef __cplusplus
extern "C" {
#endif

enum LogLevel {msg_debug, msg_verbose, msg_info, msg_warn, msg_error, msg_fatal, msg_silent};

void msg_set_loglevel(const enum LogLevel level);
void msg_set_prefix(const char prefix[]);

void msg_printf(const enum LogLevel level, char const * const fmt, ...);
void msg_abort(char const * const fmt, ...);

#ifdef __cplusplus
}
#endif
#endif
