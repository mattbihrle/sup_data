library(microeco)
library(file2meco)
 |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> |> ion and its growth rate (r).
p2 <- ggplot(data, aes(x = Abundance, y = Growth_Rate, color = Species)) +
  geom_point(size = 4) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Growth Rate (r) vs. Abundance", x = "Abundance", y = "Growth Rate (r)")

# Print the plots
print(p1)
print(p2)
|> |> |> 