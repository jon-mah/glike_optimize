# sftp_instructions

Goal: Provide some means by which the Michigan Imputation team can download via SFTP Bryan's data.

## Non-solution 1: Provide SFTP site information and access to their team

We **cannot** do this, as CARC only allows for SFTP download to local machines when we have **USC NetID, password, and duo authentication**.  
* [Source 1 (using an SFTP GUI)](https://www.carc.usc.edu/user-guides/research-data-management/transferring-data/transferring-data-gui). Quote: "Along with USC login credentials, Duo two-factor authentication (2FA) is required for CARC access."
* [Source 2 (using SFTP via command-line)](https://www.carc.usc.edu/user-guides/research-data-management/transferring-data/transferring-data-command-line.html) Quote: "To use sftp in interactive mode, from your local computer, first log in to a CARC node like hpc-transfer1 and authenticate via Duo:"

## Potential solution 1: We upload to them using SFTP

It is [recommended](https://www.carc.usc.edu/user-guides/research-data-management/transferring-data/transferring-data-command-line.html) that we use `module avail lftp` and the [`lftp`](https://lftp.yar.ru/) software to transfer files to SFTP-like servers and protocols.

Else, their server likely has some dedicated method of file upload using SFTP, which we can access (it's just that we can't have outside users access the CARC SFTP server)

## Potential solution 2: Ask ITS for an MFT

We can use USC ITS (information technology services) MFT (managed file transfer) to facilitate secure file transfers between CARC and an "external business partner", of which the Michigan Imputation server might count. **It sounds like this is mostly something that occurs on ITS/MFT's end**, however the SFTP endpoint information is as follows:
* Environment: Production
* Protocol: SFTP
* Endpoint: sftp://sftp.usc.edu:22

Server host keypoints are:
* SHA256:2DH09cQWt+wQXpIQZYbEU98QvEb6HQ6cnTuUbf0TN1Q sftp.usc.edu (RSA)
* SHA256:n0HzVaR0OkLY0ADFaOVw6O+VJmYwzEJagOm+ZRKX4tQ sftp.usc.edu (ECDSA)
* SHA256:1dv/I2wfxDA8Zp2TgX6k74UaH0XyjjrFx0S1QBEjpao sftp.usc.edu (ED25519)

Authentication methods are:
* Username + SSH Key (only RSA, ECDSA, ED25519 cryptographic algorithms)
* Username + Password
(Authentication using SSH key is strongly reccomended)

Business partner must provide:
* An SSH Public key for authentication (if using username + SSH Key)
* The egress public IP(s) of the partner’s connection to USC for the onboarding team to allow-list on USC’s perimeter firewalls

USC ITS MFT must provide:
* MFT Account username
* MFT Account password (if using username + password authentication)
* SFTP Endpoints (available in this document)
* SFTP Server Host Key Fingerprints (available in this document)