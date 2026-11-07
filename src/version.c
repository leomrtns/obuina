/* SPDX-License-Identifier: GPL-3.0-or-later */
#include "obuina.h"
const char *bt_version(void) { return BT_VERSION; }
const char *bt_status_string(crpx_status status) { return crpx_status_string(status); }
