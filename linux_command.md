## Linux Server Usage Guide & Remote Access

This document provides **basic guidance for using the shared Linux workstation** and explains how to **connect remotely via SSH**.  
The server is a **shared computational resource**, and responsible usage is essential to ensure stability, data integrity, and fair access for all users.

Before using the server, please carefully read the **caution and resource usage guidelines** below.

---

## Basic Linux Command Line

The following figure summarizes commonly used Linux command-line operations (file navigation, permissions, compression, process monitoring, etc.), which are essential for working on the server.

<p align="center">
  <img width="1000"
       alt="Basic Linux command line cheatsheet"
       src="https://github.com/user-attachments/assets/265ece83-4d22-446b-9ef1-4e7af3b6bedc" />
</p>

If you are unfamiliar with Linux, please ensure you understand:
- Directory navigation (`cd`, `ls`, `pwd`)
- File manipulation (`cp`, `mv`, `rm`, `mkdir`)
- Permissions (`chmod`, `chown`)
- Process monitoring (`top`, `htop`, `ps`)
- Disk usage (`df -h`, `du -sh`)

---

## How to Connect to the Workstation Remotely (SSH)

The workstation can be accessed remotely using **SSH**.

- **Server IP**: `141.223.143.45`
- **Access policy**:  
  Each user must use **their own account**.  
  Please request your **username and temporary password** from **Jiyun Kim**.

### SSH Connection Command

```bash
ssh ${user_name}@141.223.143.45
```
On first login, you may be asked to confirm the host key fingerprint—this is expected.  
💡 Tip: If you frequently access the server, consider using SSH keys instead of passwords for better security.  

## ⚠️ Warning ⚠️
### Important Usage Warnings (Please Read Carefully)
This server is shared among multiple users. Misuse can affect everyone.  

### File & Directory Safety
* DO NOT delete or modify other users’ files or directories  
* DO NOT run rm -rf unless you are absolutely certain of the path  
* Always double-check paths before deleting files  

### Accidental deletion on shared systems is often irreversible.

## Software & Environment Management
* DO NOT remove or modify system-wide installed packages  
* DO NOT update system libraries without permission  
* Use conda/mamba environments or local installations in your home directory when possible  

## Resource Usage Policy
* The workstation has finite shared resources:  
* CPU: up to 28 cores  
* Memory (RAM): up to 256 GB  
* Please follow these rules:  
* ❌ Do NOT use all CPU threads (-t 28, -p 28, etc.)  
* ❌ Do NOT allocate excessive memory unless explicitly approved  
* ✅ Limit thread usage (e.g., 4–8 threads unless otherwise instructed)  
* ✅ Monitor usage with top or htop  
* Excessive CPU or RAM usage can crash running jobs from other users.  

## Long-running Jobs  
* Avoid running heavy jobs directly in an interactive SSH session  
* Use nohup, tmux, or screen for long analyses  
* Clearly name output directories and logs for traceability  

## Storage Considerations  
* Regularly check your disk usage:  
```bash
du -sh ~/*
```
* Remove temporary or intermediate files when no longer needed  
* Avoid duplicating large datasets unnecessarily

## General Best Practices  
* Keep your work organized within your own home directory  
* Document commands and parameters for reproducibility  
* Log resource-heavy jobs  
* When unsure, ask before running  

* If you encounter: Login issues / Disk space problems / Performance bottlenecks / Permission errors  
Please contact Jiyun Kim before attempting system-level changes.

---

Please use the server responsibly—otherwise you might find yourself accidentally starting a PhD in our lab :) 
