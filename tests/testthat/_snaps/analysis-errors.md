# analysis and database messages keep their layout

    Code
      calculate_tfr(rates)
    Condition
      Error in `calculate_tfr()`:
      ! Age bands in `rates` must meet without gaps or overlaps.
      x The bands are [15, 20) and [22, 25).

---

    Code
      get_mortality_hmd("USA", 2000, "user", "pass")
    Condition
      Error in `get_mortality_hmd()`:
      ! Couldn't download "Mx_1x1" for "USA" from HMD.
      i Check the country code, item, credentials and connection.
      Caused by error in `HMDHFDplus::readHMDweb()`:
      ! HTTP 401

