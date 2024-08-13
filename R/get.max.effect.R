get.max.effect_tau <- function(predictResult)
{
  apply(
    predictResult, 1, function(x) 
    {
      Beta0 = x[4]
      Beta1 = x[5]
      Beta2 = x[6]
      if (Beta2 < 0 && (-Beta1 / (2 * Beta2)) >= 0 && (-Beta1 / (2 * Beta2)) <= 1)
      {
        tau = Beta0 - ((Beta1*Beta1) / (4 * Beta2));
        gmax = -Beta1 / (2 * Beta2);
      }
      else
      {
        if (Beta0 > (Beta0 + Beta1 + Beta2))
        {
          tau = Beta0;
          gmax = 0;
        }
        else
        {
          tau = (Beta0 + Beta1 + Beta2);
          gmax = 1;
        }
      }
      tau
    }
  )
}

get.max.effect_gmax <- function(predictResult)
{
  apply(
    predictResult, 1, function(x) 
    {
      Beta0 = x[4]
      Beta1 = x[5]
      Beta2 = x[6]
      if (Beta2 < 0 && (-Beta1 / (2 * Beta2)) >= 0 && (-Beta1 / (2 * Beta2)) <= 1)
      {
        tau = Beta0 - ((Beta1*Beta1) / (4 * Beta2));
        gmax = -Beta1 / (2 * Beta2);
      }
      else
      {
        if (Beta0 > (Beta0 + Beta1 + Beta2))
        {
          tau = Beta0;
          gmax = 0;
        }
        else
        {
          tau = (Beta0 + Beta1 + Beta2);
          gmax = 1;
        }
      }
      gmax
    }
  )
}

get.max.effect <- function(predictResult)
{
  gmax <- get.max.effect_gmax(predictResult)
  tau <- get.max.effect_tau(predictResult)
  data.frame(gmax,tau)
}