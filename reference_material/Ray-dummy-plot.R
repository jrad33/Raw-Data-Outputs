pdf("/Users/ray/VT_Research_Raymond_Lee/R/inflow_outflow_timeseries_leak_experiment.pdf", width = 6, height = 4)
par(oma = c(0, 2, 1, 2))
par(mar = c(2, 2, 0, 0))
par(mgp = c(3, 0.7, 0))
par(xaxs = "i", yaxs = "i")
par(xpd = FALSE)

# not counting the few days before with ET (before plastic tarp was put on)
plot(x = NULL, y = NULL, axes = FALSE, xlab = "", ylab = "", 
     xlim = c(range(inflow.agg.L.day[ , 1][start.datetime:end.datetime])), 
     ylim = c(0, 230))

par(xpd = NA)
points(inflow.agg.L.day[ , 1][start.datetime:end.datetime], inflow.agg.L.day[ , 2][start.datetime:end.datetime], pch = 19, type = "o")
points(outflow.agg.L.day[ , 1][start.datetime:end.datetime], outflow.agg.L.day[ , 2][start.datetime:end.datetime], pch = 21, bg = "lightgrey", type = "o")
par(xpd = FALSE)

axis(1, at = c(range(outflow.agg.L.day[ , 1][start.datetime:end.datetime])), lwd.ticks = 0, labels = c("",""), lwd = 2)
axis(1, at = seq(head(as.Date(outflow.agg.L.day[ , 1][start.datetime:end.datetime], origin = "1970-01-01"), n = 1), tail(as.Date(outflow.agg.L.day[ , 1][start.datetime:end.datetime], origin = "1970-01-01"), n = 1), by = "weeks"), 
     lwd = 0, lwd.ticks = 2, tcl = -0.5, 
     labels = format(seq(head(as.Date(outflow.agg.L.day[ , 1][start.datetime:end.datetime], origin = "1970-01-01"), n = 1), tail(as.Date(outflow.agg.L.day[ , 1][start.datetime:end.datetime], origin = "1970-01-01"), n = 1), by = "weeks"), "%b %d"))
axis(2, at = c(0, 225), lwd.ticks = 2, labels = c("",""), lwd = 2)
axis(2, at = seq(0, 225, by = 25), lwd = 0, lwd.ticks = 2, tcl = -0.5)

legend("topleft", c("Input (Qin)", "Total outflow (Qout)"), pch = c(19, 21), lty = c(1, 1), 
       lwd = c(2, 2), pt.bg = c(NA, "lightgrey"), inset = 0.05,
       text.width = c(5, 5), col = c("black", "black"), bty = "n")

mtext("Time ( d )", side = 1, line = 2, cex = 1.5, outer = TRUE)
mtext(bquote("Q ( " ~ L ~ d^-1 ~ ")"), side = 2, line = 0, cex = 1.5, outer = TRUE)
dev.off()