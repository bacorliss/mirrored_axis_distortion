# Relative change is a true measure of change
# RC = (Y-X)/X

# Define conversion table between rc values x<0 and x>=0
prcs <- c(.1, .25, .5, .75, 1, 1.25, 1.5, 2, 3, 4, 5, 10, 20, 30, 50, 100)
df_rc = data.frame(prc=prcs)
df_rc$nrc = mirror_rc(-df_rc$prc, forward = FALSE)



# Testing if conversion functions work
rc_test <- data.frame(rc = c(rev(df_rc$nrc), 1, df_rc$prc))
rc_test$rc_2_mrc <- mirror_rc(rc_test$rc)
rc_test$mrc_2_rc <- mirror_rc(rc_test$rc_2_mrc, forward = FALSE)


