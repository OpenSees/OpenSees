//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
#if 0
static const char*
classic_banner = R"(

         OpenSees -- Open System For Earthquake Engineering Simulation
                 Pacific Earthquake Engineering Research Center
                        Version 3.4.0 64-Bit

      (c) Copyright 1999-2022 The Regents of the University of California
                         All Rights Reserved
  (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)

)";

static const char*
peer_banner = R"(
    OpenSees -- Open System For Earthquake Engineering Simulation
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
            Pacific Earthquake Engineering Research Center
)";
#endif

static const char*
unicode_banner 
#if 1
= R"(
                 ┌─┐┌─┐┌─┐┌─┐  ┌──┌─┐┌─┐ ┌──
                 └─┘├─┘└──┘ │ ─┘  └──└───┘
 ───────────────────┘Berkeley, California ──────────────────────
                         © UC Regents
)";

#else
= R"(
                 ╔═╗╔═╗╔═╗╔═╗  ╔══╔═╗╔═╗ ╔══
                 ╚═╝╠═╝╚══╝ ║ ═╝  ╚══╚═══╝
 ═══════════════════╝Berkeley, California ══════════════════════
)";
#endif

static const char*
copyright = R"(
Copyright (c) 1999-2023 The Regents of the University of California.
All Rights Reserved.
)";

static const char*
license = R"(
Copyright @ 1999-2022 The Regents of the University of California (The
Regents). All Rights Reserved.

The Regents grants permission, without fee and without a written license
agreement, for (a) use, reproduction, modification, and distribution of this
software and its documentation by educational, research, and non-profit
entities for noncommercial purposes only; and (b) use, reproduction and
modification of this software by other entities for internal purposes only. The
above copyright notice, this paragraph and the following three paragraphs must
appear in all copies and modifications of the software and/or documentation.


Permission to incorporate this software into products for commercial
distribution may be obtained by contacting the University of California

Office of Technology Licensing
2150 Shattuck Avenue #510
Berkeley, CA 94720-1620
(510) 643-7201

This software program and documentation are copyrighted by The Regents of the
University of California. The Regents does not warrant that the operation of
the program will be uninterrupted or error-free. The end-user understands that
the program was developed for research purposes and is advised not to rely
exclusively on the program for any reason.

IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL,
INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF REGENTS HAS BEEN
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  REGENTS GRANTS NO EXPRESS OR
IMPLIED LICENSE IN ANY PATENT RIGHTS OF REGENTS BUT HAS IMPLEMENTED AN
INDIVIDUAL CONTRIBUTOR LICENSE AGREEMENT FOR THE OPENSEES PROJECT AT THE
UNIVERISTY OF CALIFORNIA, BERKELEY TO BENEFIT THE END USER.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
)";

