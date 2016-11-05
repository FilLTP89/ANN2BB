#!/bin/bash
wine exsim
wine rspectra < rspectra.dat
wine med_sp < med_sp_sa.dat
wine rspectra < rspectra_sd.dat
wine med_sp < med_sp_sd.dat
wine integration < int.dat
