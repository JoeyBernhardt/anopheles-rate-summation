


#### Plotting the posteriors from Kerri's analyses to make figures for the paper

cons_briere_bite_rate <- load("data-raw/posteriors/constant_br.c.briere_T.uniform.Rdata")
bite.rate.constant.preds <-model.out$BUGSoutput$sims.list$z.trait.mu.pred

### something like this:

predictions <- as.data.frame(cons_briere_bite_rate$BUGSoutput$sims.list$z.trait.mu.pred, col_names = Temp.xs)  ### columns are temperatures, rows are iterations
colnames(predictions) <- Temp.xs

predictions_sub <- predictions %>% 
	mutate(iteration = rownames(.)) %>% 
	select(iteration, everything()) %>% 
	sample_n(size = 1000)
